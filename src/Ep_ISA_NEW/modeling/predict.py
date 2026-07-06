import numpy as np
from loguru import logger
from Ep_ISA_NEW.utils import one_hot_encode


def _detect_input_format(model):
    """
    Detect whether model expects (N, 4, L) or (N, L, 4).
    Returns 'channels_first' or 'channels_last'.
    """
    candidates = []
    if hasattr(model, 'input_shape') and model.input_shape:
        candidates.append(model.input_shape)
    if hasattr(model, 'inputs') and model.inputs:
        for inp in model.inputs:
            candidates.append(inp.shape)
    if hasattr(model, 'layers') and model.layers:
        first = model.layers[0]
        if hasattr(first, 'input_shape') and first.input_shape:
            candidates.append(first.input_shape)

    for ishape in candidates:
        shape = tuple(s for s in (ishape if isinstance(ishape, (list, tuple)) else [ishape]))
        if len(shape) == 3:
            if shape[1] == 4:
                return 'channels_first'
            if shape[2] == 4:
                return 'channels_last'

    return 'channels_last'


def compute_predictions(model, seqs, device=None, batch_size=1024, tracks=[0]):
    """
    TF/Keras inference for DNA sequences. Auto-detects input format
    (channels_first/last) and output format (single/multi-task/dict).

    device param is accepted but ignored (TF manages GPU placement).
    """
    x_np = one_hot_encode(seqs)

    input_format = _detect_input_format(model)
    if input_format == 'channels_last':
        x_np = np.transpose(x_np, (0, 2, 1))
    elif input_format == 'channels_first':
        logger.info("Model uses channels_first (N,4,L) input")

    x_np = x_np.astype('float32')

    all_preds = []
    for i in range(0, len(seqs), batch_size):
        batch_x = x_np[i : i + batch_size]
        preds = model(batch_x, training=False)

        if isinstance(preds, dict):
            keys = sorted(preds.keys())
            preds = np.column_stack([
                p.numpy() if hasattr(p, 'numpy') else np.asarray(p)
                for p in [preds[k] for k in keys]
            ])
        elif isinstance(preds, (list, tuple)):
            preds = np.column_stack([
                p.numpy() if hasattr(p, 'numpy') else np.asarray(p)
                for p in preds
            ])
        else:
            preds = preds.numpy() if hasattr(preds, 'numpy') else np.asarray(preds)

        all_preds.append(preds)

    result = np.concatenate(all_preds, axis=0)

    if result.ndim == 1:
        result = result.reshape(-1, 1)

    if result.shape[1] == 1:
        return result

    valid_tracks = [t for t in tracks if t < result.shape[1]]
    if not valid_tracks:
        logger.warning(f"All tracks {tracks} out of range for "
                       f"model with {result.shape[1]} outputs")
        return result

    return result[:, valid_tracks]
