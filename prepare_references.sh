#!/usr/bin/env bash
set -euo pipefail

usage() {
	cat <<EOF
Usage: $0 --ref REF_NAME --dest DEST_ROOT [--sheet REFERENCE_SHEET]

Options:
	-r|--ref    Reference key to prepare (hg38|chm13|mm39|mm10)
	-d|--dest   Destination root directory where references are stored
	-s|--sheet  Path to reference sheet TSV (default: reference_sheet.tsv)
	-h|--help   Show this help and exit
EOF
}

canonicalize() {
	local p="$1"
	if command -v realpath >/dev/null 2>&1; then
		realpath "$p" 2>/dev/null && return 0
	fi
	if command -v readlink >/dev/null 2>&1; then
		readlink -f "$p" 2>/dev/null && return 0
	fi
	if command -v python3 >/dev/null 2>&1; then
		python3 -c "import os,sys;print(os.path.realpath(sys.argv[1]))" "$p" 2>/dev/null && return 0
	fi
	printf "%s" "$p"
}

download_and_prepare() {
	local ref_key="$1" dest_root="$2" sheet_path="$3"

	local url="" fasta_name="" tmp_archive=""
	case "$ref_key" in
		hg38)
			url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
			fasta_name="hg38.fasta"
			;;
		chm13)
			url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz"
			fasta_name="chm13.fasta"
			;;
		mm39)
			url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/seqs_for_alignment_pipelines/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz"
			fasta_name="mm39.fasta"
			;;
		mm10)
			url="https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/initial/mm10.fa.gz"
			fasta_name="mm10.fasta"
			;;
		*)
			echo "Unknown ref: $ref_key" >&2
			return 2
			;;
	esac

	local dest_dir="${dest_root%/}/$ref_key"
	mkdir -p "$dest_dir"

	tmp_archive="$dest_dir/$(basename "$url")"

	if command -v wget >/dev/null 2>&1; then
		wget -O "$tmp_archive" "$url"
	elif command -v curl >/dev/null 2>&1; then
		curl -L -o "$tmp_archive" "$url"
	else
		echo "Neither wget nor curl found" >&2
		return 3
	fi

	# Uncompress to fasta
	if [[ "$tmp_archive" == *.gz ]]; then
		gunzip -c "$tmp_archive" > "$dest_dir/$fasta_name"
	else
		mv "$tmp_archive" "$dest_dir/$fasta_name"
	fi

	# Index fasta
	if ! command -v samtools >/dev/null 2>&1; then
		echo "samtools not found in PATH; please install samtools to index fasta" >&2
		return 4
	fi
	samtools faidx "$dest_dir/$fasta_name"

	# Remove archive if still present
	[[ -f "$tmp_archive" ]] && rm -f "$tmp_archive"

	# Update reference sheet
	if [[ -z "$sheet_path" ]]; then
		sheet_path="reference_sheet.tsv"
	fi

	if [[ ! -f "$sheet_path" ]]; then
		printf "ref_name\tref_fasta\tref_index\n" > "$sheet_path"
	fi

	local fasta_path="$dest_dir/$fasta_name"
	local index_path="${fasta_path}.fai"

	# Use canonical (real) paths for workflow tools that dislike symlinks
	local fasta_real
	local index_real
	fasta_real=$(canonicalize "$fasta_path") || fasta_real="$fasta_path"
	index_real=$(canonicalize "$index_path") || index_real="$index_path"

	printf "%s\t%s\t%s\n" "$ref_key" "$fasta_real" "$index_real" >> "$sheet_path"

	echo "Prepared $ref_key -> $fasta_path (index: $index_path)"
}

# CLI parsing
if [[ ${#@} -eq 0 ]]; then
	usage
	exit 1
fi

REF=""
DEST=""
SHEET="reference_sheet.tsv"

while [[ $# -gt 0 ]]; do
	case "$1" in
		-r|--ref)
			REF="$2"; shift 2;;
		-d|--dest)
			DEST="$2"; shift 2;;
		-s|--sheet)
			SHEET="$2"; shift 2;;
		-h|--help)
			usage; exit 0;;
		--)
			shift; break;;
		*)
			echo "Unknown option: $1" >&2; usage; exit 1;;
	esac
done

if [[ -z "$REF" || -z "$DEST" ]]; then
	echo "--ref and --dest are required" >&2
	usage
	exit 1
fi

# Check if ref already exists in sheet
if [[ -f "$SHEET" ]] && grep -P "^$REF\t" "$SHEET"; then
	echo "Warning: $REF already exists in $SHEET" >&2
	echo "Skipping download and preparation."
	exit 0
else
	echo "Preparing reference $REF..."
fi

download_and_prepare "$REF" "$DEST" "$SHEET"


