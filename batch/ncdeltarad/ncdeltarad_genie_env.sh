# Sourced on the worker (larbatch <initsource>) by the GEN stage only.
# Points GENIE at the enhanced-radiative PDG table. The table travels with the
# job via jobsub `-f` (lands in $CONDOR_DIR_INPUT); NCDR_TABLE_PNFS is the
# fallback copied with ifdh. Later stages must NOT set GENIE_PDG_TABLE: the
# CAF stage's GENIE reweighting reads the table too.
#
# The orchestrator fills in NCDR_TABLE_PNFS and NCDR_BR_GAMMA when it stages
# this file.

NCDR_TABLE_NAME=ncdeltarad_pdg_table.txt
NCDR_TABLE_PNFS="@NCDR_TABLE_PNFS@"
NCDR_BR_GAMMA="@NCDR_BR_GAMMA@"

if [ -n "${CONDOR_DIR_INPUT:-}" ] && [ -f "${CONDOR_DIR_INPUT}/${NCDR_TABLE_NAME}" ]; then
  cp "${CONDOR_DIR_INPUT}/${NCDR_TABLE_NAME}" "./${NCDR_TABLE_NAME}"
elif command -v ifdh >/dev/null 2>&1; then
  ifdh cp "${NCDR_TABLE_PNFS}" "./${NCDR_TABLE_NAME}"
fi

if [ ! -s "./${NCDR_TABLE_NAME}" ]; then
  echo "ncdeltarad: could not obtain ${NCDR_TABLE_NAME}; refusing to generate a nominal-BR sample" >&2
  return 1
fi

# The patched Nγ lines must be present (p γ for Δ+, n γ for Δ0).
NCDR_NPATCH=$(awk -v br="${NCDR_BR_GAMMA}" '$3+0 == br+0 && $4 == 2 && $6 == 22 && ($5 == 2212 || $5 == 2112)' "./${NCDR_TABLE_NAME}" | wc -l)
if [ "${NCDR_NPATCH}" -ne 2 ]; then
  echo "ncdeltarad: ${NCDR_TABLE_NAME} does not carry BR(Nγ)=${NCDR_BR_GAMMA} for both Δ+ and Δ0" >&2
  return 1
fi

export GENIE_PDG_TABLE="${PWD}/${NCDR_TABLE_NAME}"
echo "ncdeltarad: GENIE_PDG_TABLE=${GENIE_PDG_TABLE} (BR(Nγ)=${NCDR_BR_GAMMA})"
