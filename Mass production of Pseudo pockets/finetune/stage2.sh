#!/bin/bash
# Stage 2: unfreeze the last 4 pocket-encoder transformer layers (+ pocket
# encoder's final layer norm) on top of Stage 1's result, at a 10x lower LR,
# so the pocket representation itself can adapt to RNA-binding-site geometry
# -- while everything below layer -4, and the whole mol side, stays frozen
# (bounds how far the model can drift from the pretrained DrugCLIP prior in
# one step). Restores Stage 1's best checkpoint by valid_bedroc, resets
# optimizer/scheduler state since the trainable parameter set changes.
set -euo pipefail
cd "$(dirname "$0")"

DRUGCLIP_REPO="/data1/home/sanchezxl/drugclip_project/DrugCLIP"
DATA_PATH="$(pwd)/data"
SAVE_DIR="$(pwd)/savedir_stage2"
STAGE1_CKPT="$(pwd)/savedir_stage1/checkpoint_best.pt"

if [ ! -f "$STAGE1_CKPT" ]; then
    echo "Stage 1 checkpoint not found at $STAGE1_CKPT -- run stage1.sh first" >&2
    exit 1
fi

mkdir -p "$SAVE_DIR"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS=1

unicore-train "$DATA_PATH" --user-dir "$DRUGCLIP_REPO/unimol" \
    --train-subset train --valid-subset valid \
    --num-workers 8 --ddp-backend=c10d \
    --task drugclip --loss in_batch_softmax --arch drugclip \
    --max-pocket-atoms 256 --dist-threshold 8.0 --recycling 3 \
    --restore-file "$STAGE1_CKPT" \
    --reset-optimizer --reset-lr-scheduler --reset-meters --reset-dataloader \
    --freeze-encoders --unfreeze-last-pocket-layers 4 \
    --optimizer adam --adam-betas "(0.9, 0.999)" --adam-eps 1e-8 --clip-norm 1.0 \
    --lr-scheduler polynomial_decay --lr 1e-4 --warmup-ratio 0.06 \
    --max-epoch 15 --batch-size 48 --batch-size-valid 128 \
    --update-freq 1 --seed 1 \
    --tensorboard-logdir "$SAVE_DIR/tsb" \
    --log-interval 100 --log-format simple \
    --validate-interval 1 \
    --best-checkpoint-metric valid_bedroc --patience 5 --all-gather-list-size 2048000 \
    --save-dir "$SAVE_DIR" --tmp-save-dir "$SAVE_DIR/tmp" --keep-last-epochs 3 --keep-best-checkpoints 1 \
    --find-unused-parameters \
    --maximize-best-checkpoint-metric \
    2>&1 | tee "$SAVE_DIR/train.log"
