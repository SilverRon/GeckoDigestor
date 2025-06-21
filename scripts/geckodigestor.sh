#!/bin/bash

# ?? ?? ??
SESSION="gecko"

# tmux ??? ?? ????? ??
tmux has-session -t $SESSION 2>/dev/null

if [ $? != 0 ]; then
  # ?? ?? ? ? ?? ??? ?? (AlertReceiver ??)
  tmux new-session -d -s $SESSION -n alert_receiver
  tmux send-keys -t $SESSION "cd ~/GeckoDigestor/src" C-m
  tmux send-keys -t $SESSION "conda activate geckodigestor" C-m
  tmux send-keys -t $SESSION "python AlertReceiver.py" C-m

  # ? ?? ??? ?? (GeckoDigestor ??)
  tmux new-window -t $SESSION -n digestor
  tmux send-keys -t $SESSION:1 "cd ~/GeckoDigestor/test" C-m
  tmux send-keys -t $SESSION:1 "conda activate geckodigestor" C-m
  tmux send-keys -t $SESSION:1 "python GeckoDigestor.py" C-m
else
  echo "tmux ?? '$SESSION'?(?) ?? ?????. 'tmux attach-session -t $SESSION'? ???? ?????."
fi

# ????? ??? ??
tmux attach-session -t $SESSION

