#!/usr/bin/env bash
#
# Install MToolBox: create (or update) the conda environment and add the
# mtoolbox-activate alias to a shell startup file.
set -euo pipefail

usage() {
    cat <<EOF
Usage: bash install.sh [-n env_name] [-r startup_file] [-s]

  -n  name of the conda environment (default: mtoolbox)
  -r  shell startup file where the mtoolbox-activate alias is written
      (default: ~/.bashrc)
  -s  skip creating/updating the conda environment
  -h  show this help
EOF
}

env_name=mtoolbox
rc_file="${HOME}/.bashrc"
skip_env=false
while getopts "n:r:sh" opt; do
    case $opt in
        n) env_name=$OPTARG ;;
        r) rc_file=$OPTARG ;;
        s) skip_env=true ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 1 ;;
    esac
done

# absolute path of the MToolBox folder, wherever the script is run from
mtoolbox_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v conda > /dev/null; then
    echo "conda not found: please install Anaconda or Miniconda first." >&2
    exit 1
fi

if [ "$skip_env" = false ]; then
    # mamba solves the environment much faster, if available
    solver=conda
    command -v mamba > /dev/null && solver=mamba
    if conda env list | awk '{print $1}' | grep -qx "$env_name"; then
        echo "Updating conda environment ${env_name}..."
        "$solver" env update -n "$env_name" -f "${mtoolbox_dir}/envs/mtoolbox.yaml" --prune
    else
        echo "Creating conda environment ${env_name}..."
        "$solver" env create -n "$env_name" -f "${mtoolbox_dir}/envs/mtoolbox.yaml"
    fi
fi

# (re)write the alias, replacing the one from a previous installation
alias_line="alias mtoolbox-activate='conda activate ${env_name} && export PATH=\"${mtoolbox_dir}:${mtoolbox_dir}/scripts:\$PATH\"'"
touch "$rc_file"
tmp_rc="$(mktemp)"
grep -v "^alias mtoolbox-activate=" "$rc_file" > "$tmp_rc" || true
echo "$alias_line" >> "$tmp_rc"
cat "$tmp_rc" > "$rc_file"
rm -f "$tmp_rc"

echo "
MToolBox installed.
The mtoolbox-activate command was added to ${rc_file}: open a new shell
(or run 'source ${rc_file}'), then run

    mtoolbox-activate

to activate the ${env_name} environment and add the MToolBox wrappers
(MToolBox-variant-calling, MToolBox-human-prioritization, ...) to your PATH."
