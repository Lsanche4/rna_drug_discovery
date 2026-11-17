#!/bin/bash
# Stage 1: linear-probe fine-tune of the pretrained DrugCLIP checkpoint on
# RBP pseudo-pockets. mol_model, pocket_model and mol_project are frozen
# (--freeze-encoders); only pocket_project + logit_scale train. This adapts
# the pocket embedding to the new (RNA-fragment) ligand distribution while
# anchoring it to the frozen, already-calibrated ligand embedding space --
# the safest possible fine-tune step, so it doubles as the sanity check that
# the whole data/task/loss pipeline runs end-to-end before Stage 2 touches
# any encoder weights.
set -euo pipefail
cd "$(dirname "$0")"

DRUGCLIP_REPO="/data1/home/sanchezxl/drugclip_project/DrugCLIP"
DATA_PATH="$(pwd)/data"
SAVE_DIR="$(pwd)/savedir_stage1"
PRETRAINED_CKPT="$DRUGCLIP_REPO/checkpoint_best.pt"

mkdir -p "$SAVE_DIR"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS=1

unicore-train "$DATA_PATH" --user-dir "$DRUGCLIP_REPO/unimol" \
    --train-subset train --valid-subset valid \
    --num-workers 8 --ddp-backend=c10d \
    --task drugclip --loss in_batch_softmax --arch drugclip \
    --max-pocket-atoms 256 --dist-threshold 8.0 --recycling 3 \
    --restore-file "$PRETRAINED_CKPT" \
    --reset-optimizer --reset-lr-scheduler --reset-meters --reset-dataloader \
    --freeze-encoders \
    --optimizer adam --adam-betas "(0.9, 0.999)" --adam-eps 1e-8 --clip-norm 1.0 \
    --lr-scheduler polynomial_decay --lr 1e-3 --warmup-ratio 0.06 \
    --max-epoch 30 --batch-size 48 --batch-size-valid 128 \
    --update-freq 1 --seed 1 \
    --tensorboard-logdir "$SAVE_DIR/tsb" \
    --log-interval 100 --log-format simple \
    --validate-interval 1 \
    --best-checkpoint-metric valid_bedroc --patience 5 --all-gather-list-size 2048000 \
    --save-dir "$SAVE_DIR" --tmp-save-dir "$SAVE_DIR/tmp" --keep-last-epochs 3 --keep-best-checkpoints 1 \
    --find-unused-parameters \
    --maximize-best-checkpoint-metric \
    2>&1 | tee "$SAVE_DIR/train.log"
