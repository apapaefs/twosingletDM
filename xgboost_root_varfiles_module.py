import gzip
import numpy as np
import random
import math
import csv
import json
import os
import re
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
# import xgboost and sklearn stuff:
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import pickle
from tqdm.auto import tqdm

# functions to load root varfiles from HwSim
from read_root_varfiles import *

###############
# FUNCTIONS   #
###############

# main training function (June 12, 2025: in progress!)
def train_xgboost(signal_SM_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, seed):

    # load signal and backgrounds
    # NOTE THAT: the weights will also be multiplied by the total cross section for the process!

    # load signal:
    idS=0 # id number for signal
    S, LS, wS = read_ROOT_varfile(signal_SM_file, idS, xsS)
    Sweight = Luminosity * np.sum(wS)/initial_S * sig_factors # calculate total expected number of events
    print('Signal pre-efficiency=', np.sum(wS)/initial_S/xsS)

    # initial values for arrays used in training:
    X = S
    L = LS
    W = wS

    Bweight = 0
    initial_NB = {}
    for bkg in Backgrounds:
        xsB=Backgrounds_xsec[(Energy, bkg)] # background cross sections (fb)
        B, LB, wB =  read_ROOT_varfile(Background_files[(Energy, bkg)], idB[bkg], Backgrounds_xsec[(Energy, bkg)])
        initial_NB[bkg] =  Luminosity * np.sum(wB)/initial_B[bkg] * bkg_factors # calculate total expected number of events in each background
        Bweight += initial_NB[bkg] # incremenet to total expected number of events
        print('Background pre-efficiency', bkg, np.sum(wB)/initial_B[bkg]/Backgrounds_xsec[(Energy, bkg)])
        # concatenate lists:
        X = X + B
        L = L + LB
        W = W + wB

    # convert to numpy arrays:
    X = np.array(X)
    L = np.array(L)
    W = np.array(W)

    # create testing and training samples:
    print("Splitting samples into testing and training")
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, L, W, test_size=0.5,random_state=seed)

    # train XGBoost model:
    print("Training the model")
    model = xgb.XGBClassifier()
    model.fit(X_train, y_train,sample_weight=w_train, verbose=3)
    print("Done training the model")
    return model

# apply the given model
def apply_xgboost(model, signal_SM_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, seed):

    # load signal:
    idS=0 # id number for signal
    S, LS, wS = read_ROOT_varfile(signal_SM_file, idS, xsS)
    Sweight = Luminosity * np.sum(wS)/initial_S * sig_factors # calculate total expected number of events
    print('Signal pre-efficiency=', np.sum(wS)/initial_S/xsS)

    # initial values for arrays used in training:
    X = S
    L = LS
    W = wS

    #print(model)
    Bweight = 0
    initial_NB = {}
    for bkg in Backgrounds:
        xsB=Backgrounds_xsec[(Energy, bkg)] # background cross sections (fb)
        print(Background_files[(Energy, bkg)])
        B, LB, wB =  read_ROOT_varfile(Background_files[(Energy, bkg)], idB[bkg], Backgrounds_xsec[(Energy, bkg)])
        initial_NB[bkg] =  Luminosity * np.sum(wB)/initial_B[bkg] * bkg_factors # calculate total expected number of events in each background
        Bweight += initial_NB[bkg] # incremenet to total expected number of events
        print('Background pre-efficiency', bkg, np.sum(wB)/initial_B[bkg]/Backgrounds_xsec[(Energy, bkg)])
        # concatenate lists:
        X = X + B
        L = L + LB
        W = W + wB


    # create testing and training samples:
    print("Splitting samples into testing and training")
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, L, W, test_size=0.5,random_state=seed)

    # make predictions for test data
    y_pred = model.predict(X_test)
    predictions = [round(value) for value in y_pred]

    # evaluate predictions
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))

    # Confusion matrix whose i-th row and j-th column entry indicates the number of samples with true label being i-th class and predicted label being j-th class.
    # in this case signal = 0, backgrounds = i = 1, 2,...
    # (0,0): signal-as-signal -> True positive
    # (i,0): background-as-signal (mis-id) -> False positive
    confmatrix = confusion_matrix(y_test, predictions)
    print('confusion matrix:')
    print(confmatrix)
    # signal efficiency:
    total_S = 0
    for j in range(len(Backgrounds)+1):
        total_S += confmatrix[0][j]
    eff_S = confmatrix[0][0]/total_S # signal identified as signal divided by total number of signal events
    # background effiencies:
    eff_B = {}
    for bkg in Backgrounds:
        total_B = 0
        for j in range(len(Backgrounds)+1):
            total_B += confmatrix[idB[bkg]][j]
        eff_B[bkg] = confmatrix[idB[bkg]][0]/total_B
        print(bkg, confmatrix[idB[bkg]][0], total_B)

    print('Luminosity=', Luminosity)

    # initial cross sections into final state:
    print('Initial signal cross section=', Sweight/Luminosity)
    print('Initial background cross section=', Bweight/Luminosity)
    print('-')
    # calculate "significance"
    print('Initial significance=', Sweight/np.sqrt(Bweight))
    print('-')
    # print analysis efficiencies
    print('Signal efficiency (xgboost only)=', eff_S)
    print('Signal efficiency full=', eff_S * np.sum(wS)/initial_S/xsS)
    print('Background Efficiencies (xgboost only)=', eff_B)
    print('-')
    print('Final signal cross section=', xsS *eff_S * np.sum(wS)/initial_S )
    # calculate the number of events for the background after the analysis:
    final_NB = {}
    final_NB_total = 0
    for bkg in Backgrounds:
        final_NB[bkg] = initial_NB[bkg] * eff_B[bkg]
        #print('\tNumber of events in', bkg,final_NB[bkg], 'after analysis')
        final_NB_total += final_NB[bkg]
    print('Final background cross section=', final_NB_total/Luminosity)
    print('Final significance=', Sweight*eff_S/np.sqrt(final_NB_total))
    print('-')
    # calculate 95% C.L. limit on expected number of events:
    S2sigma = np.sqrt(final_NB_total) * 2
    print('95% C.L. limit on number of signal events=', S2sigma)
    print('95% C.L. limit on signal cross section in given final state=', S2sigma/Luminosity, 'fb')


# save the model:
def save_model(model, filename):
    model.save_model(str(filename))
    #with open(filename,'wb') as f:
    #    pickle.dump(model,f)

# load the model:
def load_model(filename):
    model = xgb.XGBClassifier()
    model.load_model(filename)
    #with open(filename, 'rb') as f:
    #    model = pickle.load(f)
    return model


def _as_path_list(paths):
    if paths is None:
        return []
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(path) for path in paths]


def _expand_per_file(values, files, default):
    if values is None:
        return [default for _ in files]
    if isinstance(values, (int, float)):
        return [float(values) for _ in files]
    expanded = [default if value is None else float(value) for value in values]
    if len(expanded) == 1 and len(files) > 1:
        return expanded * len(files)
    if len(expanded) != len(files):
        raise ValueError(f"Expected {len(files)} values, got {len(expanded)}")
    return expanded


def _expand_metadata(metadata, files):
    if metadata is None:
        return [{} for _ in files]
    metadata = [dict(item or {}) for item in metadata]
    if len(metadata) == 1 and len(files) > 1:
        return metadata * len(files)
    if len(metadata) != len(files):
        raise ValueError(f"Expected {len(files)} metadata rows, got {len(metadata)}")
    return metadata


def _balanced_training_weights(labels, raw_weights):
    labels = np.asarray(labels)
    raw_weights = np.asarray(raw_weights, dtype=float)
    base = np.abs(raw_weights)
    base[base == 0.0] = 1.0

    balanced = np.zeros_like(base, dtype=float)
    for label in np.unique(labels):
        mask = labels == label
        class_sum = np.sum(base[mask])
        if class_sum > 0.0:
            balanced[mask] = base[mask] * np.sum(mask) / class_sum
        else:
            balanced[mask] = 1.0
    return balanced


def _normalisation_denominator(generated, raw_weights, normalisation_weight=None):
    if normalisation_weight is not None and normalisation_weight > 0:
        return float(normalisation_weight), "input_weight_sum"
    if generated is not None and generated > 0:
        return float(generated), "generated_events"
    fallback = float(np.sum(np.abs(raw_weights))) if np.sum(np.abs(raw_weights)) > 0 else float(len(raw_weights))
    return fallback, "loaded_weight_sum"


def _load_signal_background_group(
    files,
    label,
    xsecs_fb,
    generated_events,
    luminosity,
    max_events,
    rate_factors=None,
    normalisation_weights=None,
    sample_metadata=None,
):
    rows = []
    labels = []
    raw_weights = []
    physical_weights = []
    sources = []
    summaries = []
    if rate_factors is None:
        rate_factors = [1.0 for _ in files]
    if normalisation_weights is None:
        normalisation_weights = [None for _ in files]
    if sample_metadata is None:
        sample_metadata = [{} for _ in files]

    for path, xsec_fb, generated, rate_factor, normalisation_weight, metadata in zip(
        files, xsecs_fb, generated_events, rate_factors, normalisation_weights, sample_metadata
    ):
        features, sample_labels, weights = read_ROOT_varfile(path, label, 1.0, max_events=max_events)
        weights = np.asarray(weights, dtype=float)
        normalisation, normalisation_source = _normalisation_denominator(generated, weights, normalisation_weight)

        effective_xsec_fb = float(xsec_fb) * float(rate_factor)
        physical = luminosity * effective_xsec_fb * weights / normalisation
        preselected_events = float(np.sum(physical))
        initial_events = float(luminosity) * effective_xsec_fb
        analysis_efficiency = preselected_events / initial_events if initial_events != 0.0 else 0.0

        rows.extend(features)
        labels.extend(sample_labels)
        raw_weights.extend(weights.tolist())
        physical_weights.extend(physical.tolist())
        sources.extend([str(path)] * len(features))
        summary = dict(metadata or {})
        summary.update(
            {
                "file": str(path),
                "entries": int(len(features)),
                "sum_weight": float(np.sum(weights)),
                "xsec_fb": float(xsec_fb),
                "raw_xsec_fb": float(xsec_fb),
                "rate_factor": float(rate_factor),
                "effective_xsec_fb": float(effective_xsec_fb),
                "generated_events": None if generated is None else int(generated),
                "normalisation_weight": float(normalisation),
                "normalisation_source": normalisation_source,
                "expected_preselected_events": preselected_events,
                "analysis_efficiency": float(analysis_efficiency),
            }
        )
        summaries.append(summary)

    return rows, labels, raw_weights, physical_weights, sources, summaries


def _mc_event_count_summary(signal_summary, background_summary):
    def total_entries(rows):
        return int(sum(int(row.get("entries", 0)) for row in rows))

    def total_generated(rows):
        values = [row.get("generated_events") for row in rows]
        known = [int(value) for value in values if value is not None]
        return int(sum(known)) if len(known) == len(values) and values else None

    def total_known_generated(rows):
        return int(sum(int(row["generated_events"]) for row in rows if row.get("generated_events") is not None))

    return {
        "signal_entries_read": total_entries(signal_summary),
        "signal_generated_events": total_generated(signal_summary),
        "signal_known_generated_events": total_known_generated(signal_summary),
        "signal_files": len(signal_summary),
        "background_entries_read": total_entries(background_summary),
        "background_generated_events": total_generated(background_summary),
        "background_known_generated_events": total_known_generated(background_summary),
        "background_files": len(background_summary),
    }


def _best_significance_threshold(scores, labels, physical_weights, systematics=0.0):
    thresholds = np.linspace(0.0, 1.0, 501)
    best = {
        "threshold": 0.5,
        "signal_events": 0.0,
        "background_events": 0.0,
        "significance": 0.0,
        "signal_efficiency": 0.0,
        "background_efficiency": 0.0,
    }

    labels = np.asarray(labels)
    scores = np.asarray(scores)
    physical_weights = np.asarray(physical_weights, dtype=float)
    total_signal = np.sum(physical_weights[labels == 1])
    total_background = np.sum(physical_weights[labels == 0])

    for threshold in thresholds:
        selected = scores >= threshold
        signal = float(np.sum(physical_weights[(labels == 1) & selected]))
        background = float(np.sum(physical_weights[(labels == 0) & selected]))
        if background <= 0.0:
            continue
        denominator = math.sqrt(background + (systematics * background) ** 2)
        significance = signal / denominator
        if significance > best["significance"]:
            best = {
                "threshold": float(threshold),
                "signal_events": signal,
                "background_events": background,
                "significance": float(significance),
                "signal_efficiency": float(signal / total_signal) if total_signal > 0 else 0.0,
                "background_efficiency": float(background / total_background) if total_background > 0 else 0.0,
            }
    return best


def _write_scores_csv(path, labels, scores, physical_weights, sources):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "score", "physical_weight", "source"])
        for label, score, weight, source in zip(labels, scores, physical_weights, sources):
            writer.writerow([int(label), float(score), float(weight), source])


def _write_roc_plot(path, labels, scores, physical_weights):
    fpr, tpr, _ = roc_curve(labels, scores, sample_weight=physical_weights)
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label="XGBoost")
    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.xlabel("Background efficiency")
    plt.ylabel("Signal efficiency")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def _write_feature_importance_plot(path, model, feature_names, top_n=20):
    importances = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(importances)[-top_n:]
    plt.figure(figsize=(7, max(4, 0.25 * len(order))))
    plt.barh(np.asarray(feature_names)[order], importances[order])
    plt.xlabel("Feature importance")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def run_signal_background_analysis(
    signal_files,
    background_files,
    output_dir="xgboost_results",
    signal_xsecs_fb=None,
    background_xsecs_fb=None,
    signal_rate_factors=None,
    background_rate_factors=None,
    signal_generated_events=None,
    background_generated_events=None,
    signal_normalisation_weights=None,
    background_normalisation_weights=None,
    signal_metadata=None,
    background_metadata=None,
    luminosity=3000.0,
    test_size=0.35,
    seed=12345,
    systematics=0.0,
    max_events=None,
    model_params=None,
):
    """Train and evaluate a binary XGBoost signal-vs-background analysis."""

    signal_files = _as_path_list(signal_files)
    background_files = _as_path_list(background_files)
    if not signal_files:
        raise ValueError("At least one signal ROOT variable file is required")
    if not background_files:
        raise ValueError("At least one background ROOT variable file is required")

    signal_xsecs_fb = _expand_per_file(signal_xsecs_fb, signal_files, 1.0)
    background_xsecs_fb = _expand_per_file(background_xsecs_fb, background_files, 1.0)
    signal_rate_factors = _expand_per_file(signal_rate_factors, signal_files, 1.0)
    background_rate_factors = _expand_per_file(background_rate_factors, background_files, 1.0)
    signal_generated_events = _expand_per_file(signal_generated_events, signal_files, None)
    background_generated_events = _expand_per_file(background_generated_events, background_files, None)
    signal_normalisation_weights = _expand_per_file(signal_normalisation_weights, signal_files, None)
    background_normalisation_weights = _expand_per_file(background_normalisation_weights, background_files, None)
    signal_metadata = _expand_metadata(signal_metadata, signal_files)
    background_metadata = _expand_metadata(background_metadata, background_files)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    s_rows, s_labels, s_raw, s_phys, s_sources, signal_summary = _load_signal_background_group(
        signal_files,
        1,
        signal_xsecs_fb,
        signal_generated_events,
        luminosity,
        max_events,
        signal_rate_factors,
        signal_normalisation_weights,
        signal_metadata,
    )
    b_rows, b_labels, b_raw, b_phys, b_sources, background_summary = _load_signal_background_group(
        background_files,
        0,
        background_xsecs_fb,
        background_generated_events,
        luminosity,
        max_events,
        background_rate_factors,
        background_normalisation_weights,
        background_metadata,
    )

    X = np.asarray(s_rows + b_rows, dtype=float)
    y = np.asarray(s_labels + b_labels, dtype=int)
    raw_weights = np.asarray(s_raw + b_raw, dtype=float)
    physical_weights = np.asarray(s_phys + b_phys, dtype=float)
    sources = np.asarray(s_sources + b_sources)
    training_weights = _balanced_training_weights(y, raw_weights)

    if len(np.unique(y)) != 2:
        raise ValueError("The training sample must contain both signal and background events")

    split = train_test_split(
        X,
        y,
        raw_weights,
        physical_weights,
        training_weights,
        sources,
        test_size=test_size,
        random_state=seed,
        stratify=y,
    )
    X_train, X_test, y_train, y_test, raw_train, raw_test, phys_train, phys_test, train_w, test_w, src_train, src_test = split

    params = {
        "objective": "binary:logistic",
        "eval_metric": "logloss",
        "n_estimators": 300,
        "max_depth": 3,
        "learning_rate": 0.05,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "random_state": seed,
        "n_jobs": 1,
    }
    if model_params:
        params.update(model_params)

    model = xgb.XGBClassifier(**params)
    model.fit(X_train, y_train, sample_weight=train_w)

    scores = model.predict_proba(X_test)[:, 1]
    predictions = (scores >= 0.5).astype(int)
    accuracy = accuracy_score(y_test, predictions)
    auc_unweighted = roc_auc_score(y_test, scores)
    auc_weighted = roc_auc_score(y_test, scores, sample_weight=phys_test)
    best = _best_significance_threshold(scores, y_test, phys_test, systematics)
    best_predictions = (scores >= best["threshold"]).astype(int)
    confmatrix = confusion_matrix(y_test, best_predictions, labels=[0, 1])

    model_file = output_dir / "signal_background_xgboost.json"
    roc_file = output_dir / "roc.png"
    feature_file = output_dir / "feature_importance.png"
    metrics_file = output_dir / "metrics.json"
    scores_file = output_dir / "scores.csv"

    save_model(model, model_file)
    _write_scores_csv(scores_file, y_test, scores, phys_test, src_test)
    _write_roc_plot(roc_file, y_test, scores, phys_test)
    _write_feature_importance_plot(feature_file, model, FEATURE_NAMES)

    top_features = sorted(
        zip(FEATURE_NAMES, model.feature_importances_),
        key=lambda item: item[1],
        reverse=True,
    )[:10]
    mc_event_counts = _mc_event_count_summary(signal_summary, background_summary)

    metrics = {
        "n_events": int(len(y)),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "test_size": float(test_size),
        "seed": int(seed),
        "luminosity_fb_inverse": float(luminosity),
        "systematics": float(systematics),
        "accuracy_threshold_0p5": float(accuracy),
        "auc_unweighted": float(auc_unweighted),
        "auc_weighted": float(auc_weighted),
        "best_threshold": best,
        "confusion_matrix_at_best_threshold": confmatrix.tolist(),
        "expected_preselected_signal_events": float(np.sum(physical_weights[y == 1])),
        "expected_preselected_background_events": float(np.sum(physical_weights[y == 0])),
        "mc_event_counts": mc_event_counts,
        "signal_files": signal_summary,
        "background_files": background_summary,
        "top_features": [[name, float(value)] for name, value in top_features],
        "outputs": {
            "model": str(model_file),
            "metrics": str(metrics_file),
            "scores": str(scores_file),
            "roc": str(roc_file),
            "feature_importance": str(feature_file),
        },
    }

    with open(metrics_file, "w") as handle:
        json.dump(metrics, handle, indent=2)

    print("XGBoost signal-vs-background analysis complete")
    print("AUC (weighted) =", metrics["auc_weighted"])
    print("Best threshold =", best["threshold"])
    print("Expected S, B =", best["signal_events"], best["background_events"])
    print("S/sqrt(B + syst^2 B^2) =", best["significance"])
    print("MC event counts")
    print("  Signal entries read =", mc_event_counts["signal_entries_read"])
    if mc_event_counts["signal_generated_events"] is None:
        print("  Signal generated events = unavailable")
    else:
        print("  Signal generated events =", mc_event_counts["signal_generated_events"])
    print("  Background entries read =", mc_event_counts["background_entries_read"])
    if mc_event_counts["background_generated_events"] is None:
        print("  Background generated events = unavailable")
    else:
        print("  Background generated events =", mc_event_counts["background_generated_events"])
    print("Wrote outputs to", output_dir)

    return {"model": model, "metrics": metrics}


_C3D4_RUN_PATTERN = re.compile(
    r"run_gg_4h_([^_/]+)_"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)_"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)


def _point_metadata_from_path(path):
    for item in [Path(path).stem] + [parent.name for parent in Path(path).parents]:
        match = _C3D4_RUN_PATTERN.search(item)
        if not match:
            continue
        try:
            return {
                "run_group": match.group(1),
                "c3": float(match.group(2)),
                "d4": float(match.group(3)),
            }
        except ValueError:
            continue
    return {"run_group": "", "c3": None, "d4": None}


def score_signal_files(
    signal_files,
    model_file,
    output_dir="xgboost_signal_scores",
    threshold=0.5,
    signal_xsecs_fb=None,
    signal_rate_factors=None,
    signal_generated_events=None,
    signal_normalisation_weights=None,
    signal_metadata=None,
    luminosity=3000.0,
    max_events=None,
    write_event_scores=False,
):
    """Apply a trained binary XGBoost model to additional signal-point ROOT files."""

    signal_files = _as_path_list(signal_files)
    if not signal_files:
        raise ValueError("At least one signal ROOT variable file is required")

    signal_xsecs_fb = _expand_per_file(signal_xsecs_fb, signal_files, 1.0)
    signal_rate_factors = _expand_per_file(signal_rate_factors, signal_files, 1.0)
    signal_generated_events = _expand_per_file(signal_generated_events, signal_files, None)
    signal_normalisation_weights = _expand_per_file(signal_normalisation_weights, signal_files, None)
    signal_metadata = _expand_metadata(signal_metadata, signal_files)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model = load_model(str(model_file))
    rows = []
    all_event_rows = []

    for path, xsec_fb, rate_factor, generated, normalisation_weight, file_metadata in zip(
        signal_files,
        signal_xsecs_fb,
        signal_rate_factors,
        signal_generated_events,
        signal_normalisation_weights,
        signal_metadata,
    ):
        features, labels, weights = read_ROOT_varfile(path, 1, 1.0, max_events=max_events)
        X = np.asarray(features, dtype=float)
        raw_weights = np.asarray(weights, dtype=float)
        if X.size == 0:
            continue

        normalisation, normalisation_source = _normalisation_denominator(
            generated,
            raw_weights,
            normalisation_weight,
        )

        effective_xsec_fb = float(xsec_fb) * float(rate_factor)
        physical_weights = luminosity * effective_xsec_fb * raw_weights / normalisation
        scores = model.predict_proba(X)[:, 1]
        selected = scores >= threshold

        preselected_events = float(np.sum(physical_weights))
        selected_events = float(np.sum(physical_weights[selected]))
        selected_error = float(np.sqrt(np.sum(np.square(physical_weights[selected]))))
        initial_events = float(luminosity) * effective_xsec_fb
        analysis_efficiency = preselected_events / initial_events if initial_events != 0.0 else 0.0
        weighted_efficiency = selected_events / preselected_events if preselected_events != 0.0 else 0.0
        final_efficiency = selected_events / initial_events if initial_events != 0.0 else 0.0
        raw_weight_efficiency = (
            float(np.sum(raw_weights[selected]) / np.sum(raw_weights))
            if np.sum(raw_weights) != 0.0
            else 0.0
        )
        metadata = _point_metadata_from_path(path)

        row = dict(file_metadata or {})
        row.update(
            {
                "file": str(path),
                "run_group": metadata["run_group"],
                "c3": metadata["c3"],
                "d4": metadata["d4"],
                "entries": int(len(raw_weights)),
                "selected_entries": int(np.sum(selected)),
                "threshold": float(threshold),
                "xsec_fb": float(xsec_fb),
                "raw_xsec_fb": float(xsec_fb),
                "rate_factor": float(rate_factor),
                "effective_xsec_fb": float(effective_xsec_fb),
                "generated_events": None if generated is None else int(generated),
                "normalisation_weight": float(normalisation),
                "normalisation_source": normalisation_source,
                "expected_preselected_events": preselected_events,
                "expected_selected_events": selected_events,
                "expected_selected_error": selected_error,
                "raw_sigma_eff_fb": float(float(xsec_fb) * weighted_efficiency),
                "effective_sigma_eff_fb": float(effective_xsec_fb * weighted_efficiency),
                "analysis_efficiency": float(analysis_efficiency),
                "xgboost_efficiency": float(weighted_efficiency),
                "final_efficiency": float(final_efficiency),
                "weighted_efficiency": float(weighted_efficiency),
                "raw_weight_efficiency": float(raw_weight_efficiency),
                "mean_score": float(np.average(scores, weights=np.abs(raw_weights))),
            }
        )
        rows.append(row)

        if write_event_scores:
            for index, (score, weight, physical_weight, is_selected) in enumerate(
                zip(scores, raw_weights, physical_weights, selected)
            ):
                all_event_rows.append(
                    {
                        "file": str(path),
                        "entry": index,
                        "score": float(score),
                        "raw_weight": float(weight),
                        "physical_weight": float(physical_weight),
                        "selected": bool(is_selected),
                    }
                )

    summary_csv = output_dir / "scored_signal_points.csv"
    summary_json = output_dir / "scored_signal_points.json"
    fieldnames = [
        "file",
        "process_id",
        "description",
        "local_lhe",
        "raw_xsec_pb",
        "b_quarks",
        "c_quarks",
        "light_jets",
        "c_mistags",
        "light_mistags",
        "run_group",
        "c3",
        "d4",
        "entries",
        "selected_entries",
        "threshold",
        "xsec_fb",
        "raw_xsec_fb",
        "rate_factor",
        "effective_xsec_fb",
        "generated_events",
        "normalisation_weight",
        "normalisation_source",
        "expected_preselected_events",
        "expected_selected_events",
        "expected_selected_error",
        "raw_sigma_eff_fb",
        "effective_sigma_eff_fb",
        "analysis_efficiency",
        "xgboost_efficiency",
        "final_efficiency",
        "weighted_efficiency",
        "raw_weight_efficiency",
        "mean_score",
    ]

    with open(summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(summary_json, "w") as handle:
        json.dump(
            {
                "model_file": str(model_file),
                "threshold": float(threshold),
                "luminosity_fb_inverse": float(luminosity),
                "points": rows,
            },
            handle,
            indent=2,
        )

    event_scores = None
    if write_event_scores:
        event_scores = output_dir / "scored_signal_events.csv"
        with open(event_scores, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["file", "entry", "score", "raw_weight", "physical_weight", "selected"],
            )
            writer.writeheader()
            writer.writerows(all_event_rows)

    print("Scored", len(rows), "signal point files")
    print("Wrote", summary_csv)
    print("Wrote", summary_json)
    if event_scores is not None:
        print("Wrote", event_scores)

    return rows


def score_background_files(
    background_files,
    model_file,
    output_dir="xgboost_background_scores",
    threshold=0.5,
    background_xsecs_fb=None,
    background_rate_factors=None,
    background_generated_events=None,
    background_normalisation_weights=None,
    background_metadata=None,
    luminosity=3000.0,
    max_events=None,
):
    """Apply a trained binary XGBoost model to background ROOT files."""

    background_files = _as_path_list(background_files)
    if not background_files:
        raise ValueError("At least one background ROOT variable file is required")

    background_xsecs_fb = _expand_per_file(background_xsecs_fb, background_files, 1.0)
    background_rate_factors = _expand_per_file(background_rate_factors, background_files, 1.0)
    background_generated_events = _expand_per_file(background_generated_events, background_files, None)
    background_normalisation_weights = _expand_per_file(background_normalisation_weights, background_files, None)
    background_metadata = _expand_metadata(background_metadata, background_files)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model = load_model(str(model_file))
    rows = []

    for path, xsec_fb, rate_factor, generated, normalisation_weight, file_metadata in zip(
        background_files,
        background_xsecs_fb,
        background_rate_factors,
        background_generated_events,
        background_normalisation_weights,
        background_metadata,
    ):
        features, labels, weights = read_ROOT_varfile(path, 0, 1.0, max_events=max_events)
        X = np.asarray(features, dtype=float)
        raw_weights = np.asarray(weights, dtype=float)
        if X.size == 0:
            continue

        normalisation, normalisation_source = _normalisation_denominator(
            generated,
            raw_weights,
            normalisation_weight,
        )

        effective_xsec_fb = float(xsec_fb) * float(rate_factor)
        physical_weights = luminosity * effective_xsec_fb * raw_weights / normalisation
        scores = model.predict_proba(X)[:, 1]
        selected = scores >= threshold

        preselected_events = float(np.sum(physical_weights))
        selected_events = float(np.sum(physical_weights[selected]))
        selected_error = float(np.sqrt(np.sum(np.square(physical_weights[selected]))))
        initial_events = float(luminosity) * effective_xsec_fb
        analysis_efficiency = preselected_events / initial_events if initial_events != 0.0 else 0.0
        weighted_efficiency = selected_events / preselected_events if preselected_events != 0.0 else 0.0
        final_efficiency = selected_events / initial_events if initial_events != 0.0 else 0.0
        raw_weight_efficiency = (
            float(np.sum(raw_weights[selected]) / np.sum(raw_weights))
            if np.sum(raw_weights) != 0.0
            else 0.0
        )

        row = dict(file_metadata or {})
        row.update(
            {
                "file": str(path),
                "entries": int(len(raw_weights)),
                "selected_entries": int(np.sum(selected)),
                "threshold": float(threshold),
                "xsec_fb": float(xsec_fb),
                "raw_xsec_fb": float(xsec_fb),
                "rate_factor": float(rate_factor),
                "effective_xsec_fb": float(effective_xsec_fb),
                "generated_events": None if generated is None else int(generated),
                "normalisation_weight": float(normalisation),
                "normalisation_source": normalisation_source,
                "expected_preselected_events": preselected_events,
                "expected_selected_events": selected_events,
                "expected_selected_error": selected_error,
                "raw_sigma_eff_fb": float(float(xsec_fb) * weighted_efficiency),
                "effective_sigma_eff_fb": float(effective_xsec_fb * weighted_efficiency),
                "analysis_efficiency": float(analysis_efficiency),
                "xgboost_efficiency": float(weighted_efficiency),
                "final_efficiency": float(final_efficiency),
                "weighted_efficiency": float(weighted_efficiency),
                "raw_weight_efficiency": float(raw_weight_efficiency),
                "mean_score": float(np.average(scores, weights=np.abs(raw_weights))),
            }
        )
        rows.append(row)

    summary_csv = output_dir / "scored_background_samples.csv"
    summary_json = output_dir / "scored_background_samples.json"
    fieldnames = [
        "file",
        "process_id",
        "description",
        "local_lhe",
        "raw_xsec_pb",
        "b_quarks",
        "c_quarks",
        "light_jets",
        "c_mistags",
        "light_mistags",
        "entries",
        "selected_entries",
        "threshold",
        "xsec_fb",
        "raw_xsec_fb",
        "rate_factor",
        "effective_xsec_fb",
        "generated_events",
        "normalisation_weight",
        "normalisation_source",
        "expected_preselected_events",
        "expected_selected_events",
        "expected_selected_error",
        "raw_sigma_eff_fb",
        "effective_sigma_eff_fb",
        "analysis_efficiency",
        "xgboost_efficiency",
        "final_efficiency",
        "weighted_efficiency",
        "raw_weight_efficiency",
        "mean_score",
    ]

    with open(summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "model_file": str(model_file),
        "threshold": float(threshold),
        "luminosity_fb_inverse": float(luminosity),
        "expected_preselected_events_total": float(sum(row["expected_preselected_events"] for row in rows)),
        "expected_selected_events_total": float(sum(row["expected_selected_events"] for row in rows)),
        "effective_sigma_eff_fb_total": float(sum(row["effective_sigma_eff_fb"] for row in rows)),
        "raw_sigma_eff_fb_total": float(sum(row["raw_sigma_eff_fb"] for row in rows)),
    }
    with open(summary_json, "w") as handle:
        json.dump({"metadata": metadata, "backgrounds": rows}, handle, indent=2)

    print("Scored", len(rows), "background files")
    print("Wrote", summary_csv)
    print("Wrote", summary_json)

    return {"metadata": metadata, "backgrounds": rows}


def _finite_float(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def _json_safe_float(value):
    value = _finite_float(value)
    return value


def _plot_pdf_path(path):
    return Path(path).with_suffix(".pdf")


def _plot_pdf_output(path):
    return None if path is None else str(_plot_pdf_path(path))


def _poisson_cdf(n_observed, mean):
    n_observed = int(n_observed)
    mean = float(mean)
    if n_observed < 0:
        return 0.0
    if mean < 0.0:
        raise ValueError("Poisson mean must be non-negative")
    if mean == 0.0:
        return 1.0
    if mean > 100.0:
        logs = [-mean + k * math.log(mean) - math.lgamma(k + 1.0) for k in range(n_observed + 1)]
        max_log = max(logs)
        return float(math.exp(max_log) * sum(math.exp(value - max_log) for value in logs))

    term = math.exp(-mean)
    total = term
    for k in range(1, n_observed + 1):
        term *= mean / k
        total += term
    return float(min(max(total, 0.0), 1.0))


def _poisson_median_observed(background_events):
    background_events = float(background_events)
    if background_events <= 0.0:
        return 0
    if background_events > 100.0:
        return max(0, int(math.floor(background_events + 1.0 / 3.0 - 0.02 / background_events)))

    term = math.exp(-background_events)
    total = term
    n_observed = 0
    while total < 0.5:
        n_observed += 1
        term *= background_events / n_observed
        total += term
    return n_observed


def _poisson_signal_upper_limit(background_events, confidence_level=0.95, method="cls", observed_events=None):
    background_events = float(background_events)
    confidence_level = float(confidence_level)
    method = str(method).lower()
    if background_events < 0.0:
        raise ValueError("background_events must be non-negative")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between 0 and 1")
    if method not in {"cls", "classical"}:
        raise ValueError("method must be 'cls' or 'classical'")

    alpha = 1.0 - confidence_level
    if observed_events is None:
        observed_events = _poisson_median_observed(background_events)
    observed_events = int(observed_events)
    if observed_events < 0:
        raise ValueError("observed_events must be non-negative")

    background_cdf = _poisson_cdf(observed_events, background_events)
    if method == "cls" and background_cdf <= 0.0:
        raise ValueError("Cannot compute CLs with zero background-only CDF")

    def tail_value(signal_events):
        signal_cdf = _poisson_cdf(observed_events, background_events + signal_events)
        if method == "cls":
            return signal_cdf / background_cdf
        return signal_cdf

    low = 0.0
    high = max(1.0, math.sqrt(background_events + observed_events + 1.0) + 3.0)
    while tail_value(high) > alpha:
        high *= 2.0
        if high > 1.0e9:
            raise RuntimeError("Failed to bracket Poisson upper limit")

    for _ in range(100):
        mid = 0.5 * (low + high)
        if tail_value(mid) > alpha:
            low = mid
        else:
            high = mid

    signal_events = 0.5 * (low + high)
    return {
        "method": method,
        "confidence_level": confidence_level,
        "alpha": alpha,
        "background_events": background_events,
        "observed_events": observed_events,
        "background_cdf": background_cdf,
        "signal_events": float(signal_events),
        "tail_value": float(tail_value(signal_events)),
    }


DEFAULT_C3D4_CHEBYSHEV_TERMS = (
    [(i, 0) for i in range(0, 7)]
    + [(i, 1) for i in range(0, 5)]
    + [(i, 2) for i in range(0, 3)]
)
DEFAULT_HHHH_XSEC_SOURCE_DIR = Path("/mnt/ssd2/Projects/4H/MG5_aMC_v3_5_15/gg_4h_c3d4")
DEFAULT_HHHH_XSEC_WIDE_RUNNUM = "3"
DEFAULT_HHHH_XSEC_EXPECTED_WIDE_RUNS = 17
DEFAULT_HHHH_PERTURBATIVITY_MH = 125.0
DEFAULT_HHHH_PERTURBATIVITY_V = 246.0
DEFAULT_HHHH_PERTURBATIVITY_LEVEL = 0.5
DEFAULT_HHHH_PERTURBATIVITY_SQRTS = np.arange(200.0, 5000.0, 10.0)
ATL_PHYS_PUB_2025_003_LABEL = r"ATL-PHYS-PUB-2025-003 (no syst., $L = 3000\,\mathrm{fb}^{-1}$)"
ATL_PHYS_PUB_2025_003_SOURCE_URL = "https://cds.cern.ch/record/2924772/files/ATL-PHYS-PUB-2025-003.pdf"
ATL_PHYS_PUB_2025_003_FIGURE = "Figure 7 black no-systematics curve"
ATL_PHYS_PUB_2025_003_NO_SYST_KAPPA34 = np.array(
    [
        (11.45514, -2.49205),
        (11.27509, 2.66879),
        (11.02789, 8.50124),
        (10.62500, 18.82291),
        (10.59327, 19.51219),
        (10.07379, 30.50548),
        (9.97491, 32.62637),
        (9.50561, 41.49876),
        (9.32482, 44.90986),
        (8.87766, 52.50972),
        (8.67547, 55.83245),
        (8.16632, 63.50301),
        (8.02538, 65.50018),
        (7.37530, 73.96607),
        (7.32881, 74.51396),
        (6.72521, 80.62920),
        (6.22786, 85.50725),
        (6.07512, 86.76211),
        (5.42503, 91.26900),
        (4.77494, 95.22800),
        (4.49749, 96.50053),
        (4.12559, 97.82609),
        (3.47550, 99.24001),
        (2.82541, 99.75256),
        (2.17532, 99.32838),
        (1.52524, 98.00283),
        (1.08102, 96.50053),
        (0.87515, 95.59915),
        (0.22506, 91.76388),
        (-0.42503, 87.38070),
        (-0.66116, 85.50725),
        (-1.07512, 81.61895),
        (-1.72447, 75.18558),
        (-1.78571, 74.51396),
        (-2.37456, 67.05550),
        (-2.64020, 63.50301),
        (-3.02465, 57.97102),
        (-3.37957, 52.50972),
        (-3.67473, 47.73772),
        (-4.02671, 41.49876),
        (-4.32482, 36.17886),
        (-4.60670, 30.50548),
        (-4.97491, 23.29445),
        (-5.14168, 19.51219),
        (-5.62500, 9.17285),
        (-5.65083, 8.50124),
        (-6.07586, -2.49205),
        (-6.27509, -7.22870),
        (-6.48908, -13.48533),
        (-6.91558, -24.49629),
        (-6.92444, -24.72605),
        (-7.05800, -28.95016),
        (-7.24174, -35.48957),
        (-7.57453, -44.94521),
        (-7.61659, -46.48286),
        (-7.90806, -57.49381),
        (-8.22462, -67.49735),
        (-8.24823, -68.48710),
        (-8.49395, -79.48038),
        (-8.78542, -90.49134),
        (-8.87470, -93.92011),
        (-9.01638, -101.48463),
        (-9.23111, -112.47791),
        (-9.46798, -123.48887),
        (-9.52479, -126.45811),
        (-9.63031, -134.48215),
        (-9.75945, -145.49311),
        (-9.86201, -156.48639),
        (-9.90186, -165.57087),
        (-9.89669, -167.47968),
        (-9.85021, -174.10746),
        (-9.78896, -177.39484),
        (-9.75649, -178.49063),
        (-9.63179, -181.58360),
        (-9.52479, -183.17427),
        (-9.40452, -184.00495),
        (-9.24292, -184.04030),
        (-8.87470, -182.25521),
        (-8.61423, -180.32874),
        (-8.40319, -178.49063),
        (-8.22462, -177.09438),
        (-7.95971, -174.46094),
        (-7.57453, -169.91870),
        (-7.37234, -167.47968),
        (-6.92444, -162.05373),
        (-6.45956, -156.48639),
        (-6.27509, -154.20643),
        (-5.62500, -146.35914),
        (-5.55121, -145.49311),
        (-4.97491, -138.26441),
        (-4.63474, -134.48215),
        (-4.32482, -130.77059),
        (-3.67473, -123.98374),
        (-3.62308, -123.48887),
        (-3.02465, -117.00248),
        (-2.50590, -112.47791),
        (-2.37456, -111.18770),
        (-1.83220, -106.29198),
        (-1.72447, -105.42595),
        (-1.14448, -101.48463),
        (-1.07512, -100.91905),
        (-0.52538, -96.96006),
        (-0.42503, -96.34146),
        (0.22506, -92.94804),
        (0.87515, -90.49134),
        (0.87884, -90.49134),
        (1.52524, -88.12301),
        (2.17532, -86.67374),
        (2.82541, -86.05515),
        (3.47550, -86.19654),
        (4.12559, -87.11559),
        (4.77494, -88.86533),
        (5.18079, -90.49134),
        (5.42503, -91.23365),
        (6.07512, -93.92011),
        (6.72521, -97.49028),
        (7.28527, -101.48463),
        (7.37530, -102.01485),
        (8.02538, -106.39802),
        (8.67547, -112.03606),
        (8.72417, -112.47791),
        (9.32482, -117.32061),
        (9.94835, -123.48887),
        (9.97491, -123.71863),
        (10.62500, -129.53341),
        (11.11349, -134.48215),
        (11.27509, -136.00212),
        (11.92518, -141.90527),
        (12.33766, -145.49311),
        (12.57527, -147.45493),
        (12.99734, -150.30046),
        (13.22535, -151.29021),
        (13.64448, -152.47437),
        (13.79649, -152.15624),
        (13.87544, -151.62602),
        (14.03040, -149.69954),
        (14.17208, -146.78332),
        (14.21930, -145.49311),
        (14.29900, -142.25875),
        (14.35655, -137.11559),
        (14.37057, -134.48215),
        (14.37131, -129.92224),
        (14.33146, -123.48887),
        (14.21488, -112.47791),
        (14.06656, -101.48463),
        (13.90791, -90.49134),
        (13.87544, -88.56487),
        (13.66293, -79.48038),
        (13.41795, -68.48710),
        (13.22535, -59.11983),
        (13.17960, -57.49381),
        (12.85714, -46.48286),
        (12.58264, -35.48957),
        (12.57527, -35.20679),
        (12.43506, -30.87664),
        (12.20632, -24.49629),
        (11.92518, -14.95228),
        (11.87057, -13.48533),
        (11.45514, -2.49205),
    ],
    dtype=float,
)


def _scale_to_chebyshev(value, value_range):
    xmin, xmax = value_range
    return (2.0 * value - xmin - xmax) / (xmax - xmin)


def _chebyshev_t(order, x):
    if order == 0:
        return np.ones_like(x, dtype=float)
    if order == 1:
        return np.asarray(x, dtype=float)

    t_prev = np.ones_like(x, dtype=float)
    t_curr = np.asarray(x, dtype=float)
    for _ in range(2, order + 1):
        t_prev, t_curr = t_curr, 2.0 * x * t_curr - t_prev
    return t_curr


def _chebyshev_row(c3, d4, terms, k3_range, k4_range):
    k3 = 1.0 + float(c3)
    k4 = 1.0 + float(d4)
    x3 = _scale_to_chebyshev(k3, k3_range)
    x4 = _scale_to_chebyshev(k4, k4_range)
    return np.asarray([_chebyshev_t(i, x3) * _chebyshev_t(j, x4) for i, j in terms], dtype=float)


def _fit_c3d4_chebyshev(rows, value_key, error_key, terms, k3_range, k4_range):
    fit_points = []
    for row in rows:
        c3 = _finite_float(row.get("c3"))
        d4 = _finite_float(row.get("d4"))
        value = _finite_float(row.get(value_key))
        if c3 is None or d4 is None or value is None:
            continue
        error = _finite_float(row.get(error_key))
        if error is None or error <= 0.0:
            error = max(abs(value) * 0.05, 1.0e-30)
        fit_points.append((c3, d4, value, error))

    if len(fit_points) < len(terms):
        return {
            "status": "skipped",
            "reason": f"need at least {len(terms)} finite c3/d4 points; got {len(fit_points)}",
            "n_points": len(fit_points),
            "n_terms": len(terms),
        }

    design = np.asarray([_chebyshev_row(c3, d4, terms, k3_range, k4_range) for c3, d4, _, _ in fit_points], dtype=float)
    values = np.asarray([point[2] for point in fit_points], dtype=float)
    errors = np.asarray([point[3] for point in fit_points], dtype=float)
    fit_design = design / errors[:, None]
    fit_values = values / errors

    coeffs, _, rank, singular_values = np.linalg.lstsq(fit_design, fit_values, rcond=None)
    predictions = np.dot(design, coeffs)
    residuals = values - predictions
    dof = max(len(values) - len(coeffs), 1)
    chi2_dof = float(np.dot(residuals / errors, residuals / errors) / dof)
    condition = float(np.linalg.cond(design))

    return {
        "status": "ok",
        "n_points": len(fit_points),
        "n_terms": len(terms),
        "rank": int(rank),
        "condition": condition,
        "chi2_dof": chi2_dof,
        "terms": [[int(i), int(j)] for i, j in terms],
        "k3_range": [float(k3_range[0]), float(k3_range[1])],
        "k4_range": [float(k4_range[0]), float(k4_range[1])],
        "coefficients": [float(value) for value in coeffs],
        "singular_values": [float(value) for value in singular_values],
        "residual_rms": float(np.sqrt(np.mean(np.square(residuals)))),
        "max_abs_residual": float(np.max(np.abs(residuals))),
    }


def _evaluate_c3d4_chebyshev(c3, d4, fit):
    terms = [tuple(term) for term in fit["terms"]]
    row = _chebyshev_row(c3, d4, terms, fit["k3_range"], fit["k4_range"])
    return float(np.dot(row, np.asarray(fit["coefficients"], dtype=float)))


def _evaluate_c3d4_chebyshev_grid(fit, c3_range, d4_range, n_c3, n_d4):
    terms = [tuple(term) for term in fit["terms"]]
    coeffs = np.asarray(fit["coefficients"], dtype=float)
    c3_values = np.linspace(float(c3_range[0]), float(c3_range[1]), int(n_c3))
    d4_values = np.linspace(float(d4_range[0]), float(d4_range[1]), int(n_d4))
    c3_grid, d4_grid = np.meshgrid(c3_values, d4_values)
    k3_grid = 1.0 + c3_grid
    k4_grid = 1.0 + d4_grid
    x3 = _scale_to_chebyshev(k3_grid, fit["k3_range"])
    x4 = _scale_to_chebyshev(k4_grid, fit["k4_range"])
    max_i = max(i for i, _ in terms)
    max_j = max(j for _, j in terms)
    t3 = [_chebyshev_t(i, x3) for i in range(max_i + 1)]
    t4 = [_chebyshev_t(j, x4) for j in range(max_j + 1)]
    values = np.zeros_like(c3_grid, dtype=float)
    for coeff, (i, j) in zip(coeffs, terms):
        values += coeff * t3[i] * t4[j]
    return c3_grid, d4_grid, values


def _write_c3d4_grid_plot(
    path,
    c3_grid,
    d4_grid,
    z_grid,
    colorbar_label,
    cl_target=None,
    scatter_rows=None,
    scatter_key=None,
    contour_label=None,
    selected_label=None,
    contour_legend_label=None,
):
    finite = np.isfinite(z_grid)
    if not np.any(finite):
        return None

    fig, ax = plt.subplots(figsize=(7, 5.5))
    z_min = float(np.nanmin(z_grid[finite]))
    z_max = float(np.nanmax(z_grid[finite]))
    if z_min == z_max:
        mesh = ax.pcolormesh(c3_grid, d4_grid, z_grid, shading="auto", cmap="viridis")
        fig.colorbar(mesh, ax=ax, label=colorbar_label)
    else:
        levels = np.linspace(z_min, z_max, 24)
        contour = ax.contourf(c3_grid, d4_grid, z_grid, levels=levels, cmap="viridis")
        fig.colorbar(contour, ax=ax, label=colorbar_label)
        if cl_target is not None and z_min <= cl_target <= z_max:
            line = ax.contour(c3_grid, d4_grid, z_grid, levels=[cl_target], colors=["crimson"], linewidths=2.0)
            label = contour_label if contour_label is not None else f"S/sqrt(B) = {cl_target:g}"
            if label != "":
                ax.clabel(line, fmt={cl_target: label}, inline=True, fontsize=9)
            if contour_legend_label is not None:
                ax.plot([], [], color="crimson", linewidth=2.0, label=contour_legend_label)

    if scatter_rows is not None:
        xs = []
        ys = []
        selected_xs = []
        selected_ys = []
        for row in scatter_rows:
            c3 = _finite_float(row.get("c3"))
            d4 = _finite_float(row.get("d4"))
            if c3 is None or d4 is None:
                continue
            xs.append(c3)
            ys.append(d4)
            if scatter_key is not None:
                value = _finite_float(row.get(scatter_key))
                if value is not None and cl_target is not None and value >= cl_target:
                    selected_xs.append(c3)
                    selected_ys.append(d4)
        if xs:
            ax.scatter(xs, ys, c="white", s=18, edgecolors="black", linewidths=0.35, alpha=0.85, label="scored points")
        if selected_xs:
            ax.scatter(
                selected_xs,
                selected_ys,
                facecolors="none",
                edgecolors="crimson",
                s=80,
                linewidths=1.2,
                label=selected_label if selected_label is not None else f"scored S/sqrt(B) >= {cl_target:g}",
            )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="best")

    ax.set_xlabel("c3")
    ax.set_ylabel("d4")
    ax.set_title(colorbar_label)
    ax.grid(alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    fig.savefig(_plot_pdf_path(path))
    plt.close(fig)
    return path


def _write_c3d4_scatter_or_contour(
    path,
    rows,
    value_key,
    colorbar_label,
    cl_target=None,
    contour_label=None,
    selected_label=None,
    contour_legend_label=None,
):
    points = []
    seen = {}
    for row in rows:
        c3 = _finite_float(row.get("c3"))
        d4 = _finite_float(row.get("d4"))
        value = _finite_float(row.get(value_key))
        if c3 is None or d4 is None or value is None:
            continue
        key = (c3, d4)
        seen.setdefault(key, []).append(value)

    for (c3, d4), values in seen.items():
        points.append((c3, d4, float(np.mean(values))))

    if not points:
        return None

    c3_values = np.asarray([point[0] for point in points], dtype=float)
    d4_values = np.asarray([point[1] for point in points], dtype=float)
    z_values = np.asarray([point[2] for point in points], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    contour_written = False

    unique_c3 = np.unique(c3_values)
    unique_d4 = np.unique(d4_values)
    z_min = float(np.min(z_values))
    z_max = float(np.max(z_values))
    can_contour = (
        len(points) >= 3
        and len(unique_c3) >= 2
        and len(unique_d4) >= 2
        and z_min != z_max
    )

    if can_contour:
        try:
            levels = np.linspace(z_min, z_max, 16)
            contour = ax.tricontourf(c3_values, d4_values, z_values, levels=levels, cmap="viridis")
            fig.colorbar(contour, ax=ax, label=colorbar_label)
            contour_written = True
            if cl_target is not None and z_min <= cl_target <= z_max:
                line = ax.tricontour(
                    c3_values,
                    d4_values,
                    z_values,
                    levels=[cl_target],
                    colors=["crimson"],
                    linewidths=2.0,
                )
                label = contour_label if contour_label is not None else f"S/sqrt(B) = {cl_target:g}"
                if label != "":
                    ax.clabel(line, fmt={cl_target: label}, inline=True, fontsize=9)
                if contour_legend_label is not None:
                    ax.plot([], [], color="crimson", linewidth=2.0, label=contour_legend_label)
        except Exception:
            contour_written = False

    if not contour_written:
        scatter = ax.scatter(c3_values, d4_values, c=z_values, cmap="viridis", s=45, edgecolors="black", linewidths=0.35)
        fig.colorbar(scatter, ax=ax, label=colorbar_label)

    if cl_target is not None:
        selected = z_values >= cl_target
        if np.any(selected):
            ax.scatter(
                c3_values[selected],
                d4_values[selected],
                facecolors="none",
                edgecolors="crimson",
                s=95,
                linewidths=1.4,
                label=selected_label if selected_label is not None else f"S/sqrt(B) >= {cl_target:g}",
            )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="best")

    ax.set_xlabel("c3")
    ax.set_ylabel("d4")
    ax.set_title(colorbar_label)
    ax.grid(alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    fig.savefig(_plot_pdf_path(path))
    plt.close(fig)
    return path


def _read_hhhh_xsec_pb(lhe_file):
    with gzip.open(lhe_file, "rt", errors="ignore") as stream:
        for line in stream:
            match = re.search(r"Integrated weight \(pb\)\s*:\s*([0-9.eE+-]+)", line)
            if match is not None:
                return float(match.group(1))
    raise RuntimeError("No Integrated weight found in " + str(lhe_file))


def _read_hhhh_xsec_error_pb(proc_dir, run_name, xsec_pb):
    html_file = Path(proc_dir) / "HTML" / run_name / "results.html"
    if html_file.exists():
        with open(html_file, errors="ignore") as handle:
            html = handle.read()
        html = html.replace("&plusmn;", "+/-").replace("&#177;", "+/-").replace("\u00b1", "+/-")
        match = re.search(
            r"<b>s=\s*([0-9.eE+-]+)\s*\+/-\s*([0-9.eE+-]+)\s*\(pb\)</b>",
            html,
        )
        if match is not None:
            return float(match.group(2))
    return max(abs(float(xsec_pb)) * 0.01, 1.0e-30)


def _read_hhhh_xsec_points(source_dir):
    proc_dir = Path(source_dir)
    event_dir = proc_dir / "Events"
    if not event_dir.exists():
        raise RuntimeError("hhhh cross-section Events directory not found: " + str(event_dir))

    prefix = "run_gg_4h_"
    points = []
    counts = {}
    for run_dir in sorted(event_dir.glob(prefix + "*")):
        if not run_dir.is_dir():
            continue
        run_name = run_dir.name
        rest = run_name[len(prefix):]
        parts = rest.split("_")
        if len(parts) != 3:
            continue
        runnum, c3_text, d4_text = parts
        try:
            c3 = float(c3_text)
            d4 = float(d4_text)
        except ValueError:
            continue

        lhe_file = run_dir / "unweighted_events.lhe.gz"
        if not lhe_file.exists():
            continue
        xsec_pb = _read_hhhh_xsec_pb(lhe_file)
        xsec_error_pb = _read_hhhh_xsec_error_pb(proc_dir, run_name, xsec_pb)
        points.append(
            {
                "c3": c3,
                "d4": d4,
                "xsec_pb": xsec_pb,
                "xsec_error_pb": xsec_error_pb,
                "runnum": runnum,
                "run_name": run_name,
            }
        )
        counts[runnum] = counts.get(runnum, 0) + 1

    if not points:
        raise RuntimeError("No completed hhhh cross-section runs found in " + str(event_dir))
    points.sort(key=lambda row: (row["c3"], row["d4"], row["runnum"]))
    counts = dict(sorted(counts.items(), key=lambda item: int(item[0]) if item[0].isdigit() else item[0]))
    return points, counts


def _make_hhhh_xsec_log_levels(ratio):
    positive = ratio[np.isfinite(ratio) & (ratio > 0.0)]
    if len(positive) == 0:
        raise RuntimeError("No positive normalized hhhh cross-section values found")
    min_value = float(np.min(positive))
    max_value = float(np.max(positive))
    lo = math.floor(math.log10(min_value))
    hi = math.ceil(math.log10(max_value))
    nlevels = min(max((hi - lo) * 4 + 1, 12), 80)
    return np.logspace(lo, hi, nlevels)


def _make_hhhh_xsec_line_levels(filled_levels):
    lo = math.floor(math.log10(filled_levels[0]))
    hi = math.ceil(math.log10(filled_levels[-1]))
    levels = []
    for power in range(lo, hi + 1):
        value = 10.0 ** power
        if filled_levels[0] <= value <= filled_levels[-1]:
            levels.append(value)
    return levels


def _format_hhhh_xsec_level(value):
    if value >= 1000.0 or value < 0.01:
        return "%.0e" % value
    if value >= 10.0:
        return "%.0f" % value
    return "%.2g" % value


def _hhhh_perturbativity_partial_wave(s, c3_grid, d4_grid):
    k3_squared = (1.0 + c3_grid) ** 2
    k4 = 1.0 + d4_grid
    mh2 = DEFAULT_HHHH_PERTURBATIVITY_MH ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        prefactor = (
            3.0
            * mh2
            * np.sqrt(s ** 2 - 4.0 * mh2 * s)
            / (32.0 * np.pi * s * (s - mh2) * DEFAULT_HHHH_PERTURBATIVITY_V ** 2)
        )
        bracket = (
            -k4 * (s - mh2)
            - 3.0 * k3_squared * mh2
            + (
                6.0
                * k3_squared
                * mh2
                * (s - mh2)
                / (s - 4.0 * mh2)
                * np.log(s / mh2 - 3.0)
            )
        )
        value = np.abs(prefactor * bracket)
    return np.where(np.isfinite(value), value, 0.0)


def _hhhh_perturbativity_grid(c3_grid, d4_grid):
    max_partial_wave = np.zeros_like(c3_grid, dtype=float)
    for sqrt_s in DEFAULT_HHHH_PERTURBATIVITY_SQRTS:
        current = _hhhh_perturbativity_partial_wave(sqrt_s ** 2, c3_grid, d4_grid)
        max_partial_wave = np.maximum(max_partial_wave, current)
    return max_partial_wave


def _atlas_phys_pub_2025_003_c3d4_curve():
    curve = np.array(ATL_PHYS_PUB_2025_003_NO_SYST_KAPPA34, dtype=float, copy=True)
    curve[:, 0] -= 1.0
    curve[:, 1] -= 1.0
    return curve


def _plot_atlas_phys_pub_2025_003_curve(ax):
    curve = _atlas_phys_pub_2025_003_c3d4_curve()
    ax.plot(
        curve[:, 0],
        curve[:, 1],
        color="blue",
        linewidth=2.0,
        label=ATL_PHYS_PUB_2025_003_LABEL,
    )
    return {
        "label": ATL_PHYS_PUB_2025_003_LABEL,
        "source": ATL_PHYS_PUB_2025_003_SOURCE_URL,
        "figure": ATL_PHYS_PUB_2025_003_FIGURE,
        "coordinate_system": "digitized in kappa3,kappa4 and plotted as c3=kappa3-1, d4=kappa4-1",
        "n_points": int(len(curve)),
        "c3_min": float(np.min(curve[:, 0])),
        "c3_max": float(np.max(curve[:, 0])),
        "d4_min": float(np.min(curve[:, 1])),
        "d4_max": float(np.max(curve[:, 1])),
    }


def _write_hhhh_xsec_limit_overlay_plot(
    path,
    source_dir,
    limit_c3_grid,
    limit_d4_grid,
    limit_signal_events_grid,
    required_signal_events,
    limit_contour_label,
    limit_legend_label,
    fit_terms,
    fit_k3_range,
    fit_k4_range,
    plot_c3_range,
    plot_d4_range,
    plot_n_c3,
    plot_n_d4,
    atlas_overlay_path=None,
):
    metadata = {
        "status": "not_run",
        "source_dir": str(source_dir),
    }
    try:
        points, counts = _read_hhhh_xsec_points(source_dir)
        fit = _fit_c3d4_chebyshev(
            points,
            "xsec_pb",
            "xsec_error_pb",
            fit_terms,
            fit_k3_range,
            fit_k4_range,
        )
        metadata.update(
            {
                "status": fit.get("status", "unknown"),
                "n_points": len(points),
                "run_counts": counts,
                "chebyshev_fit": fit,
            }
        )
        if fit.get("status") != "ok":
            metadata["reason"] = fit.get("reason", "unknown reason")
            return None, metadata

        wide_count = counts.get(DEFAULT_HHHH_XSEC_WIDE_RUNNUM, 0)
        if wide_count < DEFAULT_HHHH_XSEC_EXPECTED_WIDE_RUNS:
            metadata["wide_run_warning"] = (
                "run "
                + DEFAULT_HHHH_XSEC_WIDE_RUNNUM
                + " has "
                + str(wide_count)
                + " completed points; expected "
                + str(DEFAULT_HHHH_XSEC_EXPECTED_WIDE_RUNS)
            )

        sigma_sm_pb = _evaluate_c3d4_chebyshev(0.0, 0.0, fit)
        if sigma_sm_pb <= 0.0 or not math.isfinite(sigma_sm_pb):
            metadata["status"] = "skipped"
            metadata["reason"] = "fitted hhhh SM cross section is not positive"
            return None, metadata

        c3_grid, d4_grid, xsec_grid_pb = _evaluate_c3d4_chebyshev_grid(
            fit,
            plot_c3_range,
            plot_d4_range,
            plot_n_c3,
            plot_n_d4,
        )
        ratio = xsec_grid_pb / sigma_sm_pb
        ratio_positive = np.ma.masked_where((ratio <= 0.0) | (~np.isfinite(ratio)), ratio)
        levels = _make_hhhh_xsec_log_levels(ratio)
        line_levels = _make_hhhh_xsec_line_levels(levels)

        perturbativity_grid = _hhhh_perturbativity_grid(c3_grid, d4_grid)
        perturbativity_contour_drawn = False
        perturbativity_min = float(np.nanmin(perturbativity_grid))
        perturbativity_max = float(np.nanmax(perturbativity_grid))
        if perturbativity_min <= DEFAULT_HHHH_PERTURBATIVITY_LEVEL <= perturbativity_max:
            perturbativity_contour_drawn = True

        finite_limit = np.isfinite(limit_signal_events_grid)
        contour_drawn = False
        limit_min = None
        limit_max = None
        if np.any(finite_limit):
            limit_min = float(np.nanmin(limit_signal_events_grid[finite_limit]))
            limit_max = float(np.nanmax(limit_signal_events_grid[finite_limit]))
            if limit_min != limit_max and limit_min <= required_signal_events <= limit_max:
                contour_drawn = True
            metadata["limit_signal_events_min"] = limit_min
            metadata["limit_signal_events_max"] = limit_max

        def draw_overlay(output_path, include_atlas_curve=False):
            fig, ax = plt.subplots(figsize=(8.2, 6.2), constrained_layout=True)
            contour = ax.contourf(
                c3_grid,
                d4_grid,
                ratio_positive,
                levels=levels,
                norm=colors.LogNorm(vmin=levels[0], vmax=levels[-1]),
                cmap="viridis",
                extend="both",
            )
            if np.any(np.isfinite(ratio) & (ratio <= 0.0)):
                ax.contourf(c3_grid, d4_grid, ratio <= 0.0, levels=[0.5, 1.5], colors=["0.75"], alpha=0.8)
            if line_levels:
                lines = ax.contour(
                    c3_grid,
                    d4_grid,
                    ratio_positive,
                    levels=line_levels,
                    colors="white",
                    linewidths=0.55,
                )
                ax.clabel(lines, fmt=_format_hhhh_xsec_level, inline=True, fontsize=10)

            if perturbativity_contour_drawn:
                ax.contour(
                    c3_grid,
                    d4_grid,
                    perturbativity_grid,
                    levels=[DEFAULT_HHHH_PERTURBATIVITY_LEVEL],
                    colors=["black"],
                    linestyles="--",
                    linewidths=1.7,
                )
                ax.plot(
                    [],
                    [],
                    color="black",
                    linestyle="--",
                    linewidth=1.7,
                    label=r"Perturbativity $|\mathrm{Re}\,a_0| = 0.5$",
                )

            if contour_drawn:
                limit_line = ax.contour(
                    limit_c3_grid,
                    limit_d4_grid,
                    limit_signal_events_grid,
                    levels=[required_signal_events],
                    colors=["crimson"],
                    linewidths=2.0,
                )
                if limit_contour_label != "":
                    ax.clabel(limit_line, fmt={required_signal_events: limit_contour_label}, inline=True, fontsize=9)
                if limit_legend_label is not None:
                    ax.plot([], [], color="crimson", linewidth=2.0, label=limit_legend_label)

            atlas_curve_metadata = None
            if include_atlas_curve:
                atlas_curve_metadata = _plot_atlas_phys_pub_2025_003_curve(ax)

            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(loc="best", fontsize=10)

            ax.plot([0.0], [0.0], marker="o", color="white", markeredgecolor="black", markersize=5)
            ax.set_xlim(plot_c3_range)
            ax.set_ylim(plot_d4_range)
            ax.set_xlabel(r"$c_3$", fontsize=18)
            ax.set_ylabel(r"$d_4$", fontsize=18)
            ax.set_title(r"$gg \to hhhh$ at 14 TeV: $\sigma(c_3,d_4)/\sigma(0,0)$", fontsize=20)
            ax.tick_params(axis="both", labelsize=15)

            cbar = fig.colorbar(contour, ax=ax)
            cbar.set_label(r"$\sigma(c_3,d_4)/\sigma(0,0)$", fontsize=18)
            cbar.ax.tick_params(labelsize=15)

            fig.savefig(output_path, dpi=220)
            fig.savefig(_plot_pdf_path(output_path))
            plt.close(fig)
            return atlas_curve_metadata

        draw_overlay(path)
        atlas_curve_metadata = None
        if atlas_overlay_path is not None:
            atlas_curve_metadata = draw_overlay(atlas_overlay_path, include_atlas_curve=True)

        finite_ratio = ratio[np.isfinite(ratio)]
        metadata.update(
            {
                "status": "ok",
                "output": str(path),
                "output_pdf": str(_plot_pdf_path(path)),
                "sigma_sm_pb": float(sigma_sm_pb),
                "ratio_min": float(np.nanmin(finite_ratio)) if finite_ratio.size else None,
                "ratio_max": float(np.nanmax(finite_ratio)) if finite_ratio.size else None,
                "limit_contour_drawn": contour_drawn,
                "perturbativity_level": DEFAULT_HHHH_PERTURBATIVITY_LEVEL,
                "perturbativity_min": perturbativity_min,
                "perturbativity_max": perturbativity_max,
                "perturbativity_contour_drawn": perturbativity_contour_drawn,
            }
        )
        if atlas_overlay_path is not None:
            metadata["atlas_overlay_output"] = str(atlas_overlay_path)
            metadata["atlas_overlay_output_pdf"] = str(_plot_pdf_path(atlas_overlay_path))
            metadata["atlas_reference_curve"] = atlas_curve_metadata
        return path, metadata
    except Exception as error:
        metadata["status"] = "skipped"
        metadata["reason"] = str(error)
        return None, metadata


def write_c3d4_limit_scan(
    scored_rows,
    output_dir="xgboost_c3d4_limit_scan",
    background_events=0.0,
    threshold=None,
    luminosity=3000.0,
    cl_target=2.0,
    poisson_confidence_level=0.95,
    poisson_method="cls",
    poisson_observed_events=None,
    systematics=0.0,
    model_file=None,
    metrics_file=None,
    fit_signal=True,
    fit_terms=None,
    fit_k3_range=(-29.0, 31.0),
    fit_k4_range=(-699.0, 701.0),
    plot_c3_range=(-30.0, 30.0),
    plot_d4_range=(-700.0, 700.0),
    plot_n_c3=301,
    plot_n_d4=301,
    xsec_overlay=True,
    xsec_source_dir=DEFAULT_HHHH_XSEC_SOURCE_DIR,
    rate_metadata=None,
):
    """Write c3/d4 efficiencies, sigma*eff fits, and 95% CL contour plots."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    background_events = float(background_events)
    cl_target = float(cl_target)
    poisson_confidence_level = float(poisson_confidence_level)
    systematics = float(systematics)
    if background_events < 0.0:
        raise ValueError("background_events must be non-negative")

    sqrt_background = math.sqrt(background_events) if background_events > 0.0 else 0.0
    gaussian_required_signal_events = cl_target * sqrt_background if background_events > 0.0 else None
    poisson_limit = _poisson_signal_upper_limit(
        background_events,
        confidence_level=poisson_confidence_level,
        method=poisson_method,
        observed_events=poisson_observed_events,
    )
    required_signal_events = float(poisson_limit["signal_events"])
    poisson_observed_source = "median_background" if poisson_observed_events is None else "user_specified"
    poisson_method_label = "CLs" if str(poisson_limit["method"]).lower() == "cls" else "classical"
    poisson_confidence_label = f"{100.0 * poisson_confidence_level:g}%"
    poisson_contour_label = ""
    luminosity_label = rf"$L = {luminosity:g}\,\mathrm{{fb}}^{{-1}}$"
    poisson_legend_label = f"Poisson {poisson_method_label} {poisson_confidence_label} CL, {luminosity_label}"
    poisson_selected_label = f"scored S >= S95 ({required_signal_events:.3g})"
    syst_denominator = (
        math.sqrt(background_events + (systematics * background_events) ** 2)
        if background_events > 0.0
        else 0.0
    )

    rows = []
    for row in scored_rows:
        signal_events = float(row.get("expected_selected_events", 0.0))
        preselected_events = float(row.get("expected_preselected_events", 0.0))
        selected_error = float(row.get("expected_selected_error", 0.0))
        xsec_fb = float(row.get("xsec_fb", 0.0))
        rate_factor = float(row.get("rate_factor", 1.0))
        effective_xsec_fb = float(row.get("effective_xsec_fb", xsec_fb * rate_factor))
        weighted_efficiency = (
            signal_events / preselected_events
            if preselected_events != 0.0
            else float(row.get("weighted_efficiency", 0.0))
        )
        raw_sigma_eff_fb = float(row.get("raw_sigma_eff_fb", xsec_fb * weighted_efficiency))
        effective_sigma_eff_fb = float(
            row.get(
                "effective_sigma_eff_fb",
                signal_events / float(luminosity) if luminosity > 0.0 else effective_xsec_fb * weighted_efficiency,
            )
        )
        effective_sigma_eff_error_fb = selected_error / float(luminosity) if luminosity > 0.0 else 0.0
        s_over_sqrt_b = signal_events / sqrt_background if sqrt_background > 0.0 else None
        significance_with_systematics = signal_events / syst_denominator if syst_denominator > 0.0 else None
        scale_to_95cl = required_signal_events / signal_events if signal_events > 0.0 else None
        xsec_95cl_fb = xsec_fb * scale_to_95cl if scale_to_95cl is not None else None
        gaussian_signal_scale_to_target = (
            gaussian_required_signal_events / signal_events
            if gaussian_required_signal_events is not None and signal_events > 0.0
            else None
        )
        gaussian_xsec_target_fb = (
            xsec_fb * gaussian_signal_scale_to_target
            if gaussian_signal_scale_to_target is not None
            else None
        )
        gaussian_excluded_by_s_over_sqrt_b = bool(s_over_sqrt_b is not None and s_over_sqrt_b >= cl_target)
        excluded_95cl = bool(signal_events >= required_signal_events)

        rows.append(
            {
                "file": row.get("file", ""),
                "run_group": row.get("run_group", ""),
                "c3": _json_safe_float(row.get("c3")),
                "d4": _json_safe_float(row.get("d4")),
                "entries": int(row.get("entries", 0)),
                "selected_entries": int(row.get("selected_entries", 0)),
                "threshold": float(row.get("threshold", threshold if threshold is not None else 0.0)),
                "xsec_fb": xsec_fb,
                "rate_factor": rate_factor,
                "effective_xsec_fb": effective_xsec_fb,
                "generated_events": row.get("generated_events"),
                "expected_preselected_events": preselected_events,
                "expected_selected_events": signal_events,
                "expected_selected_error": selected_error,
                "raw_sigma_eff_fb": raw_sigma_eff_fb,
                "effective_sigma_eff_fb": effective_sigma_eff_fb,
                "effective_sigma_eff_error_fb": effective_sigma_eff_error_fb,
                "weighted_efficiency": float(weighted_efficiency),
                "raw_weight_efficiency": float(row.get("raw_weight_efficiency", 0.0)),
                "mean_score": float(row.get("mean_score", 0.0)),
                "background_events": background_events,
                "limit_method": f"poisson_{poisson_limit['method']}",
                "poisson_confidence_level": poisson_confidence_level,
                "poisson_observed_events": int(poisson_limit["observed_events"]),
                "required_signal_events_95cl": required_signal_events,
                "s_over_sqrt_b": _json_safe_float(s_over_sqrt_b),
                "significance_with_systematics": _json_safe_float(significance_with_systematics),
                "signal_scale_to_95cl": _json_safe_float(scale_to_95cl),
                "xsec_95cl_fb": _json_safe_float(xsec_95cl_fb),
                "excluded_95cl": excluded_95cl,
                "gaussian_cl_target_s_over_sqrt_b": cl_target,
                "gaussian_required_signal_events": _json_safe_float(gaussian_required_signal_events),
                "gaussian_signal_scale_to_target": _json_safe_float(gaussian_signal_scale_to_target),
                "gaussian_xsec_target_fb": _json_safe_float(gaussian_xsec_target_fb),
                "gaussian_excluded_by_s_over_sqrt_b": gaussian_excluded_by_s_over_sqrt_b,
                "fitted_effective_sigma_eff_fb": None,
                "fitted_expected_selected_events": None,
                "fitted_s_over_sqrt_b": None,
                "fitted_significance_with_systematics": None,
                "fitted_excluded_95cl": None,
            }
        )

    if fit_terms is None:
        fit_terms = DEFAULT_C3D4_CHEBYSHEV_TERMS
    else:
        fit_terms = [tuple(term) for term in fit_terms]

    limit_csv = output_dir / "c3d4_limit_scan.csv"
    limit_json = output_dir / "c3d4_limit_scan.json"
    fit_json = output_dir / "c3d4_sigma_eff_chebyshev_fit.json"
    cl_plot = output_dir / "c3d4_95cl_region.png"
    cl_points_plot = output_dir / "c3d4_95cl_points.png"
    efficiency_plot = output_dir / "c3d4_efficiency.png"
    sigma_eff_plot = output_dir / "c3d4_sigma_eff_fit.png"
    hhhh_xsec_overlay_plot = output_dir / "c3d4_hhhh_xsec_with_95cl.png"
    hhhh_xsec_atlas_overlay_plot = output_dir / "c3d4_hhhh_xsec_with_95cl_atl_phys_pub_2025_003.png"

    fit = None
    fit_metadata = {"status": "disabled"}
    hhhh_xsec_overlay_metadata = {"status": "disabled" if not xsec_overlay else "not_run"}
    grid_outputs = {}
    if fit_signal:
        fit_metadata = _fit_c3d4_chebyshev(
            rows,
            "effective_sigma_eff_fb",
            "effective_sigma_eff_error_fb",
            fit_terms,
            fit_k3_range,
            fit_k4_range,
        )
        if fit_metadata.get("status") == "ok":
            fit = fit_metadata
            for row in rows:
                c3 = _finite_float(row.get("c3"))
                d4 = _finite_float(row.get("d4"))
                if c3 is None or d4 is None:
                    continue
                fitted_sigma_eff_fb = _evaluate_c3d4_chebyshev(c3, d4, fit)
                fitted_signal_events = luminosity * max(fitted_sigma_eff_fb, 0.0)
                fitted_s_over_sqrt_b = fitted_signal_events / sqrt_background if sqrt_background > 0.0 else None
                fitted_significance_with_systematics = (
                    fitted_signal_events / syst_denominator if syst_denominator > 0.0 else None
                )
                row["fitted_effective_sigma_eff_fb"] = _json_safe_float(fitted_sigma_eff_fb)
                row["fitted_expected_selected_events"] = _json_safe_float(fitted_signal_events)
                row["fitted_s_over_sqrt_b"] = _json_safe_float(fitted_s_over_sqrt_b)
                row["fitted_significance_with_systematics"] = _json_safe_float(fitted_significance_with_systematics)
                row["fitted_excluded_95cl"] = bool(fitted_signal_events >= required_signal_events)

            c3_grid, d4_grid, sigma_eff_grid = _evaluate_c3d4_chebyshev_grid(
                fit,
                plot_c3_range,
                plot_d4_range,
                plot_n_c3,
                plot_n_d4,
            )
            sigma_eff_grid_clipped = np.clip(sigma_eff_grid, 0.0, None)
            fitted_signal_events_grid = luminosity * sigma_eff_grid_clipped
            cl_plot_path = _write_c3d4_grid_plot(
                cl_plot,
                c3_grid,
                d4_grid,
                fitted_signal_events_grid,
                f"Fitted selected signal events at {luminosity:g} fb^-1",
                cl_target=required_signal_events,
                scatter_rows=rows,
                scatter_key="expected_selected_events",
                contour_label=poisson_contour_label,
                selected_label=poisson_selected_label,
                contour_legend_label=poisson_legend_label,
            )
            sigma_eff_plot_path = _write_c3d4_grid_plot(
                sigma_eff_plot,
                c3_grid,
                d4_grid,
                sigma_eff_grid_clipped,
                "Fitted signal sigma x efficiency [fb]",
            )
            hhhh_xsec_overlay_path = None
            hhhh_xsec_atlas_overlay_path = None
            if xsec_overlay:
                hhhh_xsec_overlay_path, hhhh_xsec_overlay_metadata = _write_hhhh_xsec_limit_overlay_plot(
                    hhhh_xsec_overlay_plot,
                    xsec_source_dir,
                    c3_grid,
                    d4_grid,
                    fitted_signal_events_grid,
                    required_signal_events,
                    poisson_contour_label,
                    poisson_legend_label,
                    fit_terms,
                    fit_k3_range,
                    fit_k4_range,
                    plot_c3_range,
                    plot_d4_range,
                    plot_n_c3,
                    plot_n_d4,
                    atlas_overlay_path=hhhh_xsec_atlas_overlay_plot,
                )
                hhhh_xsec_atlas_overlay_path = hhhh_xsec_overlay_metadata.get("atlas_overlay_output")
                if hhhh_xsec_overlay_path is None:
                    print("hhhh cross-section overlay skipped:", hhhh_xsec_overlay_metadata.get("reason", "unknown reason"))
            grid_outputs = {
                "cl_plot": None if cl_plot_path is None else str(cl_plot_path),
                "sigma_eff_plot": None if sigma_eff_plot_path is None else str(sigma_eff_plot_path),
                "hhhh_xsec_overlay_plot": None if hhhh_xsec_overlay_path is None else str(hhhh_xsec_overlay_path),
                "hhhh_xsec_atlas_overlay_plot": hhhh_xsec_atlas_overlay_path,
            }
        else:
            print("Chebyshev fit skipped:", fit_metadata.get("reason", "unknown reason"))

    fieldnames = [
        "file",
        "run_group",
        "c3",
        "d4",
        "entries",
        "selected_entries",
        "threshold",
        "xsec_fb",
        "rate_factor",
        "effective_xsec_fb",
        "generated_events",
        "expected_preselected_events",
        "expected_selected_events",
        "expected_selected_error",
        "raw_sigma_eff_fb",
        "effective_sigma_eff_fb",
        "effective_sigma_eff_error_fb",
        "weighted_efficiency",
        "raw_weight_efficiency",
        "mean_score",
        "background_events",
        "limit_method",
        "poisson_confidence_level",
        "poisson_observed_events",
        "required_signal_events_95cl",
        "s_over_sqrt_b",
        "significance_with_systematics",
        "signal_scale_to_95cl",
        "xsec_95cl_fb",
        "excluded_95cl",
        "gaussian_cl_target_s_over_sqrt_b",
        "gaussian_required_signal_events",
        "gaussian_signal_scale_to_target",
        "gaussian_xsec_target_fb",
        "gaussian_excluded_by_s_over_sqrt_b",
        "fitted_effective_sigma_eff_fb",
        "fitted_expected_selected_events",
        "fitted_s_over_sqrt_b",
        "fitted_significance_with_systematics",
        "fitted_excluded_95cl",
    ]
    with open(limit_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    cl_points_path = _write_c3d4_scatter_or_contour(
        cl_points_plot if fit is not None else cl_plot,
        rows,
        "expected_selected_events",
        f"Scored selected signal events at {luminosity:g} fb^-1",
        cl_target=required_signal_events,
        contour_label=poisson_contour_label,
        contour_legend_label=poisson_legend_label,
        selected_label=f"S >= S95 ({required_signal_events:.3g})",
    )
    if fit is None:
        grid_outputs["cl_plot"] = None if cl_points_path is None else str(cl_points_path)
    else:
        grid_outputs["cl_points_plot"] = None if cl_points_path is None else str(cl_points_path)

    efficiency_plot_path = _write_c3d4_scatter_or_contour(
        efficiency_plot,
        rows,
        "weighted_efficiency",
        "Signal efficiency at SM-optimized XGBoost cut",
    )

    outputs = {
        "limit_csv": str(limit_csv),
        "limit_json": str(limit_json),
        "fit_json": str(fit_json),
        "cl_plot": grid_outputs.get("cl_plot"),
        "cl_plot_pdf": _plot_pdf_output(grid_outputs.get("cl_plot")),
        "cl_points_plot": grid_outputs.get("cl_points_plot"),
        "cl_points_plot_pdf": _plot_pdf_output(grid_outputs.get("cl_points_plot")),
        "sigma_eff_plot": grid_outputs.get("sigma_eff_plot"),
        "sigma_eff_plot_pdf": _plot_pdf_output(grid_outputs.get("sigma_eff_plot")),
        "hhhh_xsec_overlay_plot": grid_outputs.get("hhhh_xsec_overlay_plot"),
        "hhhh_xsec_overlay_plot_pdf": _plot_pdf_output(grid_outputs.get("hhhh_xsec_overlay_plot")),
        "hhhh_xsec_atlas_overlay_plot": grid_outputs.get("hhhh_xsec_atlas_overlay_plot"),
        "hhhh_xsec_atlas_overlay_plot_pdf": _plot_pdf_output(grid_outputs.get("hhhh_xsec_atlas_overlay_plot")),
        "efficiency_plot": None if efficiency_plot_path is None else str(efficiency_plot_path),
        "efficiency_plot_pdf": _plot_pdf_output(efficiency_plot_path),
    }
    metadata = {
        "model_file": None if model_file is None else str(model_file),
        "metrics_file": None if metrics_file is None else str(metrics_file),
        "threshold": None if threshold is None else float(threshold),
        "luminosity_fb_inverse": float(luminosity),
        "background_events": background_events,
        "limit_method": f"poisson_{poisson_limit['method']}",
        "poisson_confidence_level": poisson_confidence_level,
        "poisson_observed_source": poisson_observed_source,
        "poisson_limit": poisson_limit,
        "required_signal_events_95cl": required_signal_events,
        "sigma_eff_target_95cl_fb": required_signal_events / float(luminosity) if luminosity > 0.0 else None,
        "cl_target_s_over_sqrt_b": cl_target,
        "gaussian_cl_target_s_over_sqrt_b": cl_target,
        "gaussian_required_signal_events": _json_safe_float(gaussian_required_signal_events),
        "background_effective_sigma_eff_fb": background_events / float(luminosity) if luminosity > 0.0 else None,
        "systematics": systematics,
        "n_points": len(rows),
        "n_excluded_95cl": sum(1 for row in rows if row["excluded_95cl"]),
        "n_fitted_excluded_95cl": sum(1 for row in rows if row["fitted_excluded_95cl"]),
        "rate_metadata": rate_metadata or {},
        "chebyshev_fit": fit_metadata,
        "hhhh_xsec_overlay": hhhh_xsec_overlay_metadata,
        "outputs": outputs,
    }
    with open(fit_json, "w") as handle:
        json.dump(fit_metadata, handle, indent=2)

    with open(limit_json, "w") as handle:
        json.dump({"metadata": metadata, "points": rows}, handle, indent=2)

    print("Wrote c3/d4 limit scan", limit_csv)
    print("Wrote c3/d4 limit metadata", limit_json)
    if outputs["fit_json"] is not None:
        print("Wrote", outputs["fit_json"])
    if outputs["cl_plot"] is not None:
        print("Wrote", outputs["cl_plot"])
    if outputs["cl_plot_pdf"] is not None:
        print("Wrote", outputs["cl_plot_pdf"])
    if outputs["cl_points_plot"] is not None:
        print("Wrote", outputs["cl_points_plot"])
    if outputs["cl_points_plot_pdf"] is not None:
        print("Wrote", outputs["cl_points_plot_pdf"])
    if outputs["sigma_eff_plot"] is not None:
        print("Wrote", outputs["sigma_eff_plot"])
    if outputs["sigma_eff_plot_pdf"] is not None:
        print("Wrote", outputs["sigma_eff_plot_pdf"])
    if outputs["hhhh_xsec_overlay_plot"] is not None:
        print("Wrote", outputs["hhhh_xsec_overlay_plot"])
    if outputs["hhhh_xsec_overlay_plot_pdf"] is not None:
        print("Wrote", outputs["hhhh_xsec_overlay_plot_pdf"])
    if outputs["hhhh_xsec_atlas_overlay_plot"] is not None:
        print("Wrote", outputs["hhhh_xsec_atlas_overlay_plot"])
    if outputs["hhhh_xsec_atlas_overlay_plot_pdf"] is not None:
        print("Wrote", outputs["hhhh_xsec_atlas_overlay_plot_pdf"])
    if efficiency_plot_path is not None:
        print("Wrote", efficiency_plot_path)
    if outputs["efficiency_plot_pdf"] is not None:
        print("Wrote", outputs["efficiency_plot_pdf"])
    print(f"Poisson {poisson_confidence_label} CL target summary")
    print("  method =", poisson_method_label)
    print("  observed n =", poisson_limit["observed_events"], f"({poisson_observed_source})")
    print("  B =", background_events, "expected events")
    print("  S95 =", required_signal_events, "expected signal events")
    if luminosity > 0.0:
        print("  sigma*eff target =", required_signal_events / float(luminosity), "fb at L =", luminosity, "fb^-1")
    if gaussian_required_signal_events is not None:
        print("Gaussian S/sqrt(B) reference")
        print("  S for S/sqrt(B) =", cl_target, "is", gaussian_required_signal_events, "expected events")
    else:
        print("  S/sqrt(B) is undefined because B <= 0")

    return {"metadata": metadata, "points": rows, "outputs": outputs}

# apply the given model
def apply_xgboost_write(modelfile, signal_file, Backgrounds, Background_files, Backgrounds_xsec, xsS, initial_S, sig_factors, initial_B, idB, bkg_factors, Luminosity, Energy, seed, smeartag):
    print('loading', modelfile)
    model = xgb.XGBClassifier()
    model.load_model(modelfile)
    print('model loaded')

    # load signal:
    idS=0 # id number for signal
    S, LS, wS = read_ROOT_varfile(signal_file, idS, xsS)
    Sweight = Luminosity * np.sum(wS)/initial_S * sig_factors # calculate total expected number of events
    #print('Signal pre-efficiency=', np.sum(wS)/initial_S/xsS)

    # initial values for arrays used in training:
    X = S
    L = LS
    W = wS

    #print(model)
    Bweight = 0
    initial_NB = {}
    preeff_B = {}
    for bkg in Backgrounds:
        xsB=Backgrounds_xsec[(Energy, bkg)] # background cross sections (fb)
        B, LB, wB =  read_ROOT_varfile(Background_files[(Energy, bkg)], idB[bkg], Backgrounds_xsec[(Energy, bkg)])
        preeff_B[bkg] = np.sum(wB)/initial_B[bkg]/Backgrounds_xsec[(Energy, bkg)]
        initial_NB[bkg] =  Luminosity * np.sum(wB)/initial_B[bkg] * bkg_factors # calculate total expected number of events in each background
        Bweight += initial_NB[bkg] # incremenet to total expected number of events
        #print('Background pre-efficiency', bkg, np.sum(wB)/initial_B[bkg]/Backgrounds_xsec[(Energy, bkg)])
        # concatenate lists:
        X = X + B
        L = L + LB
        W = W + wB


    # create testing and training samples:
    #print("Splitting samples into testing and training")
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, L, W, test_size=0.99,random_state=seed)

    # make predictions for test data
    y_pred = model.predict(X_test)
    predictions = [round(value) for value in y_pred]

    # evaluate predictions
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))

    # Confusion matrix whose i-th row and j-th column entry indicates the number of samples with true label being i-th class and predicted label being j-th class.
    # in this case signal = 0, backgrounds = i = 1, 2,...
    # (0,0): signal-as-signal -> True positive
    # (i,0): background-as-signal (mis-id) -> False positive
    confmatrix = confusion_matrix(y_test, predictions)
    #print('confusion matrix:')
    #print(confmatrix)
    # signal efficiency:
    total_S = 0
    for j in range(len(Backgrounds)+1):
        total_S += confmatrix[0][j]
    eff_S = confmatrix[0][0]/total_S # signal identified as signal divided by total number of signal events
    # background effiencies:
    eff_B = {}
    for bkg in Backgrounds:
        total_B = 0
        for j in range(len(Backgrounds)+1):
            total_B += confmatrix[idB[bkg]][j]
        eff_B[bkg] = confmatrix[idB[bkg]][0]/total_B
        print(bkg, confmatrix[idB[bkg]][0], total_B)

    #print('Luminosity=', Luminosity)

    # initial cross sections into final state:
    #print('Initial signal cross section=', Sweight/Luminosity)
    #print('Initial background cross section=', Bweight/Luminosity)
    #print('-')
    # calculate "significance"
    #print('Initial significance=', Sweight/np.sqrt(Bweight))
    #print('-')
    # print analysis efficiencies
    #print('Signal efficiency=', eff_S)
    #print('Background Efficiencies=', eff_B)
    #print('-')
    #print('Final signal cross section=', Sweight/Luminosity*eff_S)
    # calculate the number of events for the background after the analysis:
    final_NB = {}
    final_NB_total = 0
    for bkg in Backgrounds:
        final_NB[bkg] = initial_NB[bkg] * eff_B[bkg]
        #print('\tNumber of events in', bkg,final_NB[bkg], 'after analysis')
        final_NB_total += final_NB[bkg]
    #print('Final background cross section=', final_NB_total/Luminosity)
    #print('Final significance=', Sweight*eff_S/np.sqrt(final_NB_total))
    #print('-')
    # calculate 95% C.L. limit on expected number of events:
    S2sigma = np.sqrt(final_NB_total) * 2
    #print('95% C.L. limit on number of signal events=', S2sigma)
    #print('95% C.L. limit on signal cross section in given final state=', S2sigma/Luminosity, 'fb')

    # open files and write all the efficiencies calculated:
    # total efficiency (including what is called "pre-efficiency"

    total_eff_S = eff_S*np.sum(wS)/initial_S/xsS
    filestream = open(signal_file.replace('_var.smear' + smeartag + '.root',  smeartag + '.XGBOOST.dat') ,'w')
    filestream.write(str(total_eff_S))
    filestream.close()
    for bkg in Backgrounds:
        filestream = open(Background_files[(Energy, bkg)].replace('_var.smear' + smeartag + '.root',  smeartag + '.XGBOOST.dat') ,'w')
        filestream.write(str(preeff_B[bkg]*eff_B[bkg]))
        filestream.close()




#########################
# Testing starts here   #
#########################

# train the model:
#trained_model, Sweight, Bweight, X_test, y_test = train_xgboost()
# save the model:
#save_model(trained_model, 'trained_model.pkl')
# laod the model:
#trained_model_test = load_model('trained_model.pkl')
# apply the model:
#apply_xgboost(trained_model_test, Sweight, Bweight, X_test, y_test)
