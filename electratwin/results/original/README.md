# Original execution failure / 原始执行失败

The unchanged extracted program exited with code 1 because this installed NumPy version has no `np.trapz`. See `execution.json`, `stdout.log` and `stderr.log`. It produced no completed scientific payload. The SyntaxWarning about `\p` in the PDE docstring is separately preserved.

Only an explicit one-line NumPy API compatibility copy was executed subsequently. Its results are stored in the sibling `compatibility` directory and retain the original scientific defects. No real instrument was connected.
