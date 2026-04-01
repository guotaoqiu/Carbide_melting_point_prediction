"""
Melting point prediction using CrabNet (local, no internet needed).

CrabNet (Compositionally-Restricted Attention-Based Network) predicts material
properties from chemical formula only. It runs locally on CPU or GPU.

Workflow:
1. Train CrabNet on a melting point dataset (formula + melting point CSV)
2. Save the trained model
3. Predict melting points for screened compounds

CrabNet does NOT ship with a pre-trained melting point model, so you must
first train it on a melting point dataset. Training data can come from:
- The curated dataset included here (data/melting_points_train.csv)
- Matminer / AFLOW databases
- Literature compilations (e.g., MAPP paper's ~10k compound dataset)

Usage:
    # Train a new model on melting point data
    python predict_melting_point_crabnet.py train --data data/melting_points_train.csv

    # Predict melting points for screened compounds
    python predict_melting_point_crabnet.py predict --input screening_all.csv --output with_mp.csv

    # Train and predict in one go
    python predict_melting_point_crabnet.py train-predict --data data/melting_points_train.csv \
        --input screening_all.csv --output with_mp.csv

References:
    - Wang et al., npj Comput. Mater. 7, 77 (2021). DOI: 10.1038/s41524-021-00545-1
    - GitHub: https://github.com/anthony-wang/CrabNet
    - Docs: https://crabnet.readthedocs.io/

Requirements:
    pip install crabnet pandas
"""

import argparse
import os

import pandas as pd


# ── Default paths ────────────────────────────────────────────────────────────

DEFAULT_MODEL_DIR = os.path.join(os.path.dirname(__file__), "models", "crabnet_melting_point")
DEFAULT_TRAIN_DATA = os.path.join(os.path.dirname(__file__), "data", "melting_points_train.csv")


def train_crabnet(
    train_csv: str,
    model_dir: str = DEFAULT_MODEL_DIR,
    epochs: int = 300,
    batch_size: int = 64,
    val_fraction: float = 0.2,
    verbose: bool = True,
):
    """
    Train CrabNet on a melting point dataset.

    The training CSV must have at least two columns:
        - formula: chemical formula string (e.g., "TiC", "LaB2C2")
        - target: melting point in Kelvin

    Args:
        train_csv: Path to training data CSV
        model_dir: Directory to save trained model
        epochs: Number of training epochs (default 300)
        batch_size: Training batch size (default 64)
        val_fraction: Fraction of data for validation (default 0.2)
        verbose: Print training progress
    """
    from crabnet.crabnet_ import CrabNet

    # Load training data
    df = pd.read_csv(train_csv)

    # Validate columns
    if "formula" not in df.columns:
        raise ValueError(f"Training CSV must have a 'formula' column. Found: {list(df.columns)}")
    if "target" not in df.columns:
        # Try common alternatives
        for alt in ["melting_point_K", "melting_point_C", "melting_temperature", "Tm"]:
            if alt in df.columns:
                df = df.rename(columns={alt: "target"})
                # Convert Celsius to Kelvin if needed
                if "C" in alt and df["target"].mean() < 4000:
                    df["target"] = df["target"] + 273.15
                break
        else:
            raise ValueError(
                f"Training CSV must have a 'target' column (melting point in K). "
                f"Found: {list(df.columns)}"
            )

    # Keep only formula + target
    df = df[["formula", "target"]].dropna()
    print(f"Training data: {len(df)} compounds")
    print(f"  Melting point range: {df['target'].min():.0f} - {df['target'].max():.0f} K")
    print(f"  Mean: {df['target'].mean():.0f} K, Std: {df['target'].std():.0f} K")

    # Initialize CrabNet
    os.makedirs(model_dir, exist_ok=True)

    cb = CrabNet(
        mat_prop="melting_point",
        learningrate=1e-3,
        batch_size=batch_size,
        epochs=epochs,
        out_dims=3,
        d_model=512,
        N=3,
        heads=4,
    )

    # Train
    print(f"\nTraining CrabNet ({epochs} epochs, batch_size={batch_size})...")
    cb.fit(train_df=df, val_df=None, test_size=val_fraction)

    # Save model
    cb.save_network(model_dir)
    print(f"Model saved to {model_dir}")

    # Report validation metrics
    val_pred = cb.predict(cb.val_df)
    if hasattr(val_pred, '__len__') and len(val_pred) > 0:
        from sklearn.metrics import mean_absolute_error, r2_score
        val_true = cb.val_df["target"].values
        if isinstance(val_pred, tuple):
            val_pred = val_pred[0]
        mae = mean_absolute_error(val_true, val_pred)
        r2 = r2_score(val_true, val_pred)
        print(f"\nValidation results:")
        print(f"  MAE:  {mae:.1f} K ({mae:.1f} C)")
        print(f"  R2:   {r2:.4f}")

    return cb


def load_crabnet(model_dir: str = DEFAULT_MODEL_DIR):
    """Load a previously trained CrabNet model."""
    from crabnet.crabnet_ import CrabNet

    cb = CrabNet(mat_prop="melting_point")
    cb.load_network(model_dir)
    print(f"Loaded CrabNet model from {model_dir}")
    return cb


def predict_melting_points(
    cb,
    input_csv: str,
    output_csv: str = "compounds_with_crabnet_mp.csv",
) -> pd.DataFrame:
    """
    Predict melting points for screened compounds using trained CrabNet.

    Args:
        cb: Trained CrabNet model instance
        input_csv: CSV from screening (must have 'formula' column)
        output_csv: Output CSV path

    Returns:
        DataFrame with added melting_point_C, mp_std_error_C, mp_source columns
    """
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} compounds from {input_csv}")

    if "formula" not in df.columns:
        raise ValueError(f"Input CSV must have a 'formula' column. Found: {list(df.columns)}")

    # Prepare prediction DataFrame (CrabNet expects 'formula' and 'target' columns)
    pred_df = df[["formula"]].copy()
    pred_df["target"] = 0  # dummy target for prediction

    # Predict
    print("Predicting melting points with CrabNet...")
    result = cb.predict(pred_df, return_uncertainty=True)

    # Parse results
    if isinstance(result, tuple) and len(result) >= 2:
        predictions_K, uncertainties = result[0], result[1]
    else:
        predictions_K = result
        uncertainties = [None] * len(df)

    # Convert K to C
    predictions_C = [round(p - 273.15) if p is not None else None for p in predictions_K]
    uncertainties_C = [round(u) if u is not None else None for u in uncertainties]

    df["melting_point_C"] = predictions_C
    df["mp_std_error_C"] = uncertainties_C
    df["mp_source"] = "crabnet"

    # Save
    df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    # Summary
    valid = df["melting_point_C"].dropna()
    print(f"\nPrediction summary:")
    print(f"  Predicted: {len(valid)}/{len(df)} compounds")
    print(f"  Range: {valid.min():.0f} - {valid.max():.0f} C")
    print(f"  Mean:  {valid.mean():.0f} C")
    if uncertainties_C[0] is not None:
        avg_unc = pd.Series(uncertainties_C).dropna().mean()
        print(f"  Avg uncertainty: {avg_unc:.0f} C")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Melting point prediction using CrabNet (local, no internet needed)"
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Train command
    train_parser = subparsers.add_parser("train", help="Train CrabNet on melting point data")
    train_parser.add_argument("--data", required=True, help="Training data CSV (formula + target columns)")
    train_parser.add_argument("--model-dir", default=DEFAULT_MODEL_DIR, help="Directory to save model")
    train_parser.add_argument("--epochs", type=int, default=300, help="Training epochs (default: 300)")
    train_parser.add_argument("--batch-size", type=int, default=64, help="Batch size (default: 64)")
    train_parser.add_argument("--val-fraction", type=float, default=0.2, help="Validation fraction (default: 0.2)")

    # Predict command
    predict_parser = subparsers.add_parser("predict", help="Predict melting points using trained model")
    predict_parser.add_argument("--input", required=True, help="Input CSV from screening")
    predict_parser.add_argument("--output", default="compounds_with_crabnet_mp.csv", help="Output CSV")
    predict_parser.add_argument("--model-dir", default=DEFAULT_MODEL_DIR, help="Directory with trained model")

    # Train + predict command
    tp_parser = subparsers.add_parser("train-predict", help="Train and predict in one step")
    tp_parser.add_argument("--data", required=True, help="Training data CSV")
    tp_parser.add_argument("--input", required=True, help="Input CSV from screening")
    tp_parser.add_argument("--output", default="compounds_with_crabnet_mp.csv", help="Output CSV")
    tp_parser.add_argument("--model-dir", default=DEFAULT_MODEL_DIR, help="Directory to save model")
    tp_parser.add_argument("--epochs", type=int, default=300, help="Training epochs")
    tp_parser.add_argument("--batch-size", type=int, default=64, help="Batch size")
    tp_parser.add_argument("--val-fraction", type=float, default=0.2, help="Validation fraction")

    args = parser.parse_args()

    if args.command == "train":
        train_crabnet(
            train_csv=args.data,
            model_dir=args.model_dir,
            epochs=args.epochs,
            batch_size=args.batch_size,
            val_fraction=args.val_fraction,
        )

    elif args.command == "predict":
        cb = load_crabnet(args.model_dir)
        predict_melting_points(cb, input_csv=args.input, output_csv=args.output)

    elif args.command == "train-predict":
        cb = train_crabnet(
            train_csv=args.data,
            model_dir=args.model_dir,
            epochs=args.epochs,
            batch_size=args.batch_size,
            val_fraction=args.val_fraction,
        )
        predict_melting_points(cb, input_csv=args.input, output_csv=args.output)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
