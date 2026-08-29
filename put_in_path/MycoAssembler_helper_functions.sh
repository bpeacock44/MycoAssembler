check_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -s "$file" ]]; then
        echo "Error: $desc ($file) not found or empty. Exiting."
        exit 1
    fi
}

validate_metadata() {
  echo "Validating metadata file."
  local line_num=0
  local err=0

  while IFS= read -r line; do
    line_num=$((line_num+1))

    # Remove leading/trailing whitespace
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"

    # Skip empty or comment lines
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    #echo "Validating line $line_num: $line"

    # split line into fields
    IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER GENOME_SIZE BUSCO_TYPE BUSCO_LIN <<< "$line"

    # Count columns
    actual_cols=$(awk -F'\t' '{print NF}' <<< "$line")
    expected_cols=11
    if (( actual_cols < expected_cols )); then
        echo "Error (line $line_num): Expected $expected_cols columns, found $actual_cols (DES=${DES:-NA})"
        err=1
        continue
    fi

    # Normalize case
    ASSEMBLER=$(tr '[:upper:]' '[:lower:]' <<< "$ASSEMBLER")
    BUSCO_TYPE=$(tr '[:upper:]' '[:lower:]' <<< "$BUSCO_TYPE")

    # ---- required fields (which can be NA) ----
    for field in LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER GENOME_SIZE BUSCO_TYPE BUSCO_LIN; do
      if [[ -z "${!field}" ]]; then
        echo "Error (line $line_num): Missing value for $field (DES=$DES)"
        err=1
      fi
    done

    # Always required (NO NA ALLOWED)
    for field in DES SPECIES ASSEMBLER BUSCO_TYPE BUSCO_LIN; do
      [[ -z "${!field}" || "${!field}" == "NA" ]] && {
        echo "Error (line $line_num): Missing value for $field (DES=$DES)"
        err=1
      }
    done

    # ---- assembler validity ----
    case "$ASSEMBLER" in
      flye|canu|spades) ;;
      *)
        echo "Error (line $line_num): Invalid assembler '$ASSEMBLER' (DES=$DES)"
        err=1
        ;;
    esac

    if [[ "$ASSEMBLER" == "flye" && "$LR" == "NA" ]]; then
      echo "Error (line $line_num): Flye requires long reads (DES=$DES)"
      err=1
    fi

    if [[ "$ASSEMBLER" == "canu" && "$LR" == "NA" ]]; then
      echo "Error (line $line_num): Canu requires long reads (DES=$DES)"
      err=1
    fi

    if [[ "$ASSEMBLER" == "spades" && "$R1" == "NA" ]]; then
      echo "Error (line $line_num): In this pipeline, SPAdes is for short reads (DES=$DES)"
      err=1
    fi

    # ---- genome size checks ----
    if [[ "$ASSEMBLER" == "canu" ]]; then
      if [[ "$GENOME_SIZE" == "NA" ]]; then
        echo "Error (line $line_num): Canu requires genome size (DES=$DES)"
        err=1
      elif ! [[ "$GENOME_SIZE" =~ ^[0-9]+[kKmMgG]$ ]]; then
        echo "Error (line $line_num): Invalid genome size '$GENOME_SIZE' (DES=$DES)"
        err=1
      fi
    fi

    # ---- BUSCO type ----
    case "$BUSCO_TYPE" in
      protein|genome|both) ;;
      *)
        echo "Error (line $line_num): Invalid BUSCO_TYPE '$BUSCO_TYPE' (DES=$DES)"
        err=1
        ;;
    esac

    # ---- input file existence ----
    [[ "$LR" != "NA" && ! -s "$LR_DIR/$LR" ]] && {
      echo "Error (line $line_num): Missing long-read file '$LR' (DES=$DES)"
      err=1
    }

    [[ "$R1" != "NA" && ! -s "$SR_DIR/$R1" ]] && {
      echo "Error (line $line_num): Missing R1 '$R1' (DES=$DES)"
      err=1
    }

    [[ "$R2" != "NA" && ! -s "$SR_DIR/$R2" ]] && {
      echo "Error (line $line_num): Missing R2 '$R2' (DES=$DES)"
      err=1
    }

    [[ "$RNA1" != "NA" && ! -s "$RNA_DIR/$RNA1" ]] && {
      echo "Error (line $line_num): Missing RNA1 '$RNA1' (DES=$DES)"
      err=1
    }

    [[ "$RNA2" != "NA" && ! -s "$RNA_DIR/$RNA2" ]] && {
      echo "Error (line $line_num): Missing RNA2 '$RNA2' (DES=$DES)"
      err=1
    }

  done < <(printf "%s\n" "$INPUT_CLEAN")

  if [[ "$err" -ne 0 ]]; then
    echo
    echo "Metadata validation failed. Fix errors above and rerun."
    exit 1
  fi
}

preprocess_reads() {
    local BASE="$1"
    local R1="$2"
    local R2="$3"
    local RNA1="$4"
    local RNA2="$5"

    local ts log_temp log_file
    ts=$(date +%Y%m%d_%H%M%S)
    log_temp=$(mktemp "$LOGDIR/part1.${BASE}.XXXXXX.log.tmp")
    log_file="$LOGDIR/part1.${BASE}.${ts}.log"

    {
        echo "--- Preprocessing reads for $BASE ---"
        echo "Started at: $(date)"
        echo

        # short read processing
        if [[ "$R1" == "NA" || "$R2" == "NA" ]]; then
            echo "--- Skipping short-read preprocessing ---"
            echo "Reason: R1 or R2 marked as NA in metadata"
            echo
        else
            local INFILE1="$SR_DIR/$R1"
            local INFILE2="$SR_DIR/$R2"

            local LEFTTRIM="$WORKDIR/${BASE}_1P.fastq.gz"
            local RIGHTTRIM="$WORKDIR/${BASE}_2P.fastq.gz"
            local LEFT="$WORKDIR/${BASE}_filtered_1.fastq.gz"
            local RIGHT="$WORKDIR/${BASE}_filtered_2.fastq.gz"

            check_file "$INFILE1" "Input left reads"
            check_file "$INFILE2" "Input right reads"

            echo "--- Step 1: Trimming short reads ---"
            if [[ ! -s "$LEFTTRIM" || ! -s "$RIGHTTRIM" ]]; then
                AAFTF trim \
                    --method bbduk \
                    --left "$INFILE1" \
                    --right "$INFILE2" \
                    -c "$THREADS" \
                    -o "$WORKDIR/$BASE"
            else
                echo "Trimmed reads already exist; skipping."
            fi

            check_file "$LEFTTRIM" "Trimmed left reads"
            check_file "$RIGHTTRIM" "Trimmed right reads"

            echo
            echo "--- Step 2: Filtering short reads ---"
            if [[ ! -s "$LEFT" || ! -s "$RIGHT" ]]; then
                AAFTF filter \
                    -c "$THREADS" \
                    --left "$LEFTTRIM" \
                    --right "$RIGHTTRIM" \
                    --aligner bbduk \
                    -o "$WORKDIR/$BASE"
            else
                echo "Filtered reads already exist; skipping."
            fi

            check_file "$LEFT" "Filtered left reads"
            check_file "$RIGHT" "Filtered right reads"
            echo
        fi

        # RNA read processing
        if [[ "$RNA1" == "NA" && "$RNA2" == "NA" ]]; then
            echo "--- Skipping RNA read preprocessing ---"
            echo "Reason: RNA1 and RNA2 both marked as NA in metadata"
            echo
        else
            echo "--- Preprocessing RNA reads ---"
        
            # Paired-end RNA
            if [[ "$RNA1" != "NA" && "$RNA2" != "NA" ]]; then
                local RNA_IN1="$SR_DIR/$RNA1"
                local RNA_IN2="$SR_DIR/$RNA2"
                local RNA_OUT1="$WORKDIR/${BASE}_RNA_1.fastq.gz"
                local RNA_OUT2="$WORKDIR/${BASE}_RNA_2.fastq.gz"
        
                check_file "$RNA_IN1" "Input RNA read 1"
                check_file "$RNA_IN2" "Input RNA read 2"
        
                echo "Detected paired-end RNA reads"
                if [[ ! -s "$RNA_OUT1" || ! -s "$RNA_OUT2" ]]; then
                    fastp \
                        -i "$RNA_IN1" \
                        -I "$RNA_IN2" \
                        -o "$RNA_OUT1" \
                        -O "$RNA_OUT2" \
                        -w "$THREADS" \
                        --detect_adapter_for_pe
                else
                    echo "Processed RNA reads already exist; skipping."
                fi
        
                check_file "$RNA_OUT1" "Filtered RNA read 1"
                check_file "$RNA_OUT2" "Filtered RNA read 2"
        
            # Single-end RNA
            else
                local RNA_IN RNA_OUT
        
                if [[ "$RNA1" != "NA" ]]; then
                    RNA_IN="$SR_DIR/$RNA1"
                    RNA_OUT="$WORKDIR/${BASE}_RNA_SE.fastq.gz"
                else
                    RNA_IN="$SR_DIR/$RNA2"
                    RNA_OUT="$WORKDIR/${BASE}_RNA_SE.fastq.gz"
                fi
        
                check_file "$RNA_IN" "Input RNA single-end read"
        
                echo "Detected single-end RNA read"
                if [[ ! -s "$RNA_OUT" ]]; then
                    fastp \
                        -i "$RNA_IN" \
                        -o "$RNA_OUT" \
                        -w "$THREADS"
                else
                    echo "Processed RNA read already exists; skipping."
                fi
        
                check_file "$RNA_OUT" "Filtered RNA single-end read"
            fi
        
            echo
        fi

        # no RNA or short reads
        if [[ ("$R1" == "NA" || "$R2" == "NA") && ("$RNA1" == "NA" || "$RNA2" == "NA") ]]; then
            echo "--- No read data available ---"
            echo "Both short-read and RNA-read preprocessing were skipped."
            echo
        fi

        echo "--- Finished preprocessing for $BASE ---"
        echo "Completed at: $(date)"

    } > "$log_temp" 2>&1

    # Strip ANSI color codes and finalize log
    sed -E 's/\x1B\[[0-9;]*[mK]//g' "$log_temp" > "$log_file"
    rm -f "$log_temp"
}


assemble_genome() {
    local BASE="$1"         # sanitized DES
    local ASSEMBLER="$2"    # flye, canu, or spades
    local LR="$3"           # long read file (or NA)
    local R1="$4"           # short read R1
    local R2="$5"           # short read R2
    local GENOME_SIZE="$6"  # e.g., 40m, 2.5g

    local log_temp log_file ASMFILE INPUT_DIR timestamp
    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part2.${BASE}.${timestamp}.log.tmp"
    log_file="$LOGDIR/part2.${BASE}.${timestamp}.log"

    # Determine input type
    if [[ "$ASSEMBLER" == "spades" ]]; then
        INPUT_DIR="$SR_DIR"
        ASMFILE="$ASM_DIR/$BASE.spades.fasta"
        # Skip if no short reads
        if [[ "$R1" == "NA" || "$R2" == "NA" ]]; then
            echo "--- Skipping SPAdes assembly for $BASE (no short reads listed)" | tee "$log_file"
            return 0
        fi
    else
        INPUT_DIR="$LR_DIR"
        if [[ "$ASSEMBLER" == "flye" ]]; then
            ASMFILE="$ASM_DIR/$BASE.flyeout/assembly.fasta"
        else
            ASMFILE="$ASM_DIR/$BASE.canuout/$BASE.contigs.fasta"
        fi
        # Skip if no long reads
        if [[ "$LR" == "NA" ]]; then
            echo "--- Skipping $ASSEMBLER assembly for $BASE (no long reads listed)" | tee "$log_file"
            return 0
        fi
    fi

    {
        echo "--- Assembly for $BASE using $ASSEMBLER ---"
        echo "Started at: $(date)"
        echo

        # Check input files
        [[ "$LR" != "NA" ]] && check_file "$INPUT_DIR/$LR" "Input long reads"
        [[ "$R1" != "NA" ]] && check_file "$INPUT_DIR/$R1" "Input left reads"
        [[ "$R2" != "NA" ]] && check_file "$INPUT_DIR/$R2" "Input right reads"

        # Run assembler if output doesn't exist
        if [[ ! -s "$ASMFILE" ]]; then
            case "$ASSEMBLER" in
                flye)
                    flye --nano-raw "$INPUT_DIR/$LR" --threads "$THREADS" -o "$ASM_DIR/$BASE.flyeout"
                    ;;
                canu)
                    canu -p "$BASE" -d "$ASM_DIR/$BASE.canuout" -nanopore "$INPUT_DIR/$LR" \
                        useGrid=false maxThreads="$THREADS" maxMemory="$MAX_MEMORY" \
                        stopOnLowCoverage=0 corOutCoverage=100 genomeSize="$GENOME_SIZE"
                    ;;
                spades)
                    AAFTF assemble --isolate -c "$THREADS" \
                        --left "$INPUT_DIR/$R1" --right "$INPUT_DIR/$R2" \
                        --memory "$MAX_MEMORY" -o "$ASMFILE" -w "$ASM_DIR/$BASE.spadesout"
                    rm -rf "$ASM_DIR/$BASE.spadesout/K??" "$ASM_DIR/$BASE.spadesout/tmp"
                    ;;
            esac
            check_file "$ASMFILE" "$ASSEMBLER assembly"
        else
            echo "Assembly already exists: $ASMFILE"
        fi

        echo
        echo "Completed at: $(date)"
    } > "$log_temp" 2>&1

    # Clean up log: strip color codes
    sed -E 's/\x1B\[[0-9;]*[mK]//g' "$log_temp" > "$log_file"
    rm -f "$log_temp"
}

polish_and_assess() {
    local BASE="$1"        # sanitized DES
    local ASSEMBLER="$2"   # flye, canu, spades
    local LR="$3"          # long read file (or NA)
    local R1="$4"          # short read R1
    local R2="$5"          # short read R2

    local LONG_FLAG=0 SHORT_FLAG=0
    [[ "$LR" != "NA" ]] && LONG_FLAG=1
    [[ "$R1" != "NA" ]] && SHORT_FLAG=1

    local ASMFILE POLISHDIR VECCLEAN MEDAKA_OUT MEDAKA_ASM MEDAKA_MODEL \
          ONT ONT_PAF POLCA SORTED STATS LEFT RIGHT log_temp log_file timestamp

    BASE="${BASE//[^a-zA-Z0-9._-]/_}"
    POLISHDIR="$ASM_DIR/$BASE.$ASSEMBLER.polishout"
    mkdir -vp "$POLISHDIR"

    if [[ "$ASSEMBLER" == "flye" ]]; then
        ASMFILE="$ASM_DIR/$BASE.flyeout/assembly.fasta"
    elif [[ "$ASSEMBLER" == "canu" ]]; then
        ASMFILE="$ASM_DIR/$BASE.canuout/$BASE.contigs.fasta"
    else
        ASMFILE="$ASM_DIR/$BASE.spades.fasta"
    fi

    VECCLEAN="$POLISHDIR/$BASE.vecscreen.fasta"
    MEDAKA_OUT="$POLISHDIR/$BASE.medakaout"
    [[ $LONG_FLAG -eq 1 ]] && mkdir -p "$MEDAKA_OUT"
    MEDAKA_ASM="$MEDAKA_OUT/consensus.fasta"
    MEDAKA_MODEL="r1041_e82_400bps_hac_v5.0.0"
    ONT="$LR_DIR/$LR"
    ONT_PAF="$MEDAKA_OUT/$BASE.paf"
    POLCA="$POLISHDIR/$BASE.polca.fasta"
    SORTED="$POLISHDIR/$BASE.sorted.fasta"
    STATS="$POLISHDIR/$BASE.sorted.stats.txt"
    LEFT="$WORKDIR/${BASE}_filtered_1.fastq.gz"
    RIGHT="$WORKDIR/${BASE}_filtered_2.fastq.gz"

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part3.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part3.$BASE.$timestamp.log"

    {
        echo "--- Post-assembly processing for $BASE ($ASSEMBLER) ---"
        echo "Threads: $THREADS | Max memory: $MAX_MEMORY GB"
        echo "Started at: $(date)"
        echo

        # Check input assembly exists
        check_file "$ASMFILE" "Assembly"

        # Step 1: Vecscreen
        if [[ ! -s "$VECCLEAN" ]]; then
            for i in $(seq 1 3); do
                echo "Running Vecscreen attempt $i/3..."
                if AAFTF vecscreen -i "$ASMFILE" -c "$THREADS" -o "$VECCLEAN"; then
                    echo "Vecscreen completed."
                    break
                else
                    echo "Vecscreen failed (attempt $i)."
                    [[ $i -lt 3 ]] && sleep 30 || { echo "Vecscreen failed 3x. Exiting."; exit 1; }
                fi
            done
        else
            echo "$VECCLEAN exists; skipping vecscreen."
        fi
        check_file "$VECCLEAN" "Vecscreen output"

        # Step 2: Medaka long-read polishing
        if (( LONG_FLAG )); then
            if [[ ! -s "$MEDAKA_ASM" ]]; then
                minimap2 -t "$THREADS" -x map-ont "$VECCLEAN" "$ONT" > "$ONT_PAF"
                medaka_consensus -i "$ONT" -d "$VECCLEAN" -o "$MEDAKA_OUT" \
                                 -t "$THREADS" -m "$MEDAKA_MODEL"
            else
                echo "Medaka output exists; skipping."
            fi
            check_file "$MEDAKA_ASM" "Medaka consensus output"
        else
            echo "No long-read data; skipping Medaka polishing."
        fi

        # Step 3: Polca short-read polishing
        if (( SHORT_FLAG )); then
            [[ $LONG_FLAG -eq 0 ]] && MEDAKA_ASM="$ASMFILE"
            if [[ ! -s "$POLCA" ]]; then
                AAFTF polish --method polca -i "$MEDAKA_ASM" -o "$POLCA" \
                             -c "$THREADS" --left "$LEFT" --right "$RIGHT" --memory "$MAX_MEMORY"
            else
                echo "$POLCA exists; skipping Polca."
            fi
            check_file "$POLCA" "Polca-polished assembly"
        else
            echo "No short-read data; skipping Polca polishing."
        fi

        # Step 4: Sort assembly
        if [[ -s "$POLCA" ]]; then
            SORT_INPUT="$POLCA"
        elif [[ -s "$MEDAKA_ASM" ]]; then
            SORT_INPUT="$MEDAKA_ASM"
        else
            SORT_INPUT=""
            echo "No polished assembly available to sort; skipping."
        fi

        if [[ -n "$SORT_INPUT" && ! -s "$SORTED" ]]; then
            echo "Sorting assembly..."
            AAFTF sort -i "$SORT_INPUT" -o "$SORTED"
        else
            echo "$SORTED exists; skipping sort."
        fi
        check_file "$SORTED" "Sorted assembly"

        # Step 5: Assess assembly
        if [[ ! -s "$STATS" ]]; then
            AAFTF assess -i "$SORTED" -r "$STATS"
        else
            echo "$STATS exists; skipping assessment."
        fi

        echo "--- Finished post-assembly for $BASE ---"
    } > "$log_temp" 2>&1

    sed -E 's/\x1B\[[0-9;]*[mK]//g' "$log_temp" > "$log_file"
    rm -f "$log_temp"
}

annotate_genome() {
    local BASE="$1"        # sanitized DES
    local ASSEMBLER="$2"   # flye, canu, spades
    local RNA1="$3"        # RNA-seq read 1 (or NA)
    local RNA2="$4"        # RNA-seq read 2 (or NA)
    local SPECIES="$5"     # species name

    local RNAseq_FLAG=0
    [[ "$RNA1" != "NA" ]] && RNAseq_FLAG=1

    BASE="${BASE//[^a-zA-Z0-9._-]/_}"
    local POLISHDIR SORTED ANNO_SAMPLE_DIR CLEAN SORT2 MASK PREDICT PREDICT_STATS ANNO_OUT

    POLISHDIR="$ASM_DIR/$BASE.$ASSEMBLER.polishout"
    SORTED="$POLISHDIR/$BASE.sorted.fasta"

    ANNO_SAMPLE_DIR="$ANNO_DIR/$BASE.$ASSEMBLER"
    mkdir -vp "$ANNO_SAMPLE_DIR"
    cd "$ANNO_SAMPLE_DIR" || exit 1

    CLEAN="$ANNO_SAMPLE_DIR/$BASE.funannotate01.fasta"
    SORT2="$ANNO_SAMPLE_DIR/$BASE.funannotate02.fasta"
    MASK="$ANNO_SAMPLE_DIR/$BASE.funannotate03.fasta"
    PREDICT="$ANNO_SAMPLE_DIR/$BASE.funannotate04"
    UPDATE="$PREDICT/update_results/${SPECIES}_${BASE}.stats.json"
    PREDICT_STATS="$PREDICT/predict_results/${SPECIES}_${BASE}.stats.json"
    ANNO_OUT="$PREDICT/annotate_results/${SPECIES}_${BASE}.stats.json"

    local timestamp log_temp log_file
    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part4.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part4.$BASE.$timestamp.log"

    {
        echo "--- Annotation for $BASE ($ASSEMBLER) ---"
        echo "Threads: $THREADS | Max memory: $MAX_MEMORY GB"
        echo "Started at: $(date)"
        
        check_file "$SORTED" "Sorted assembly input"

        echo "--- Step 1: funannotate clean ---"
        [[ ! -s "$CLEAN" ]] && funannotate clean -i "$SORTED" -o "$CLEAN" || echo "$CLEAN exists; skipping clean."

        check_file "$CLEAN" "Cleaned scaffolds"

        echo "--- Step 2: funannotate sort ---"
        [[ ! -s "$SORT2" ]] && funannotate sort -i "$CLEAN" -o "$SORT2" --minlen 50 --base scaffold || echo "$SORT2 exists; skipping sort."
        check_file "$SORT2" "Sorted scaffolds"

        echo "--- Step 3: funannotate mask ---"
        [[ ! -s "$MASK" ]] && funannotate mask -i "$SORT2" -o "$MASK" --cpus "$THREADS" || echo "$MASK exists; skipping mask."
        check_file "$MASK" "Masked scaffolds"

        if (( RNAseq_FLAG )); then
            echo "--- Step 4: funannotate train ---"
            if [[ ! -s "$PREDICT/training/funannotate_train.pasa.gff3" ]]; then
                if [[ "$RNA2" != "NA" ]]; then
                    funannotate train -i "$MASK" -o "$PREDICT" \
                        --left "${RNA_DIR}/${RNA1}" --right "${RNA_DIR}/${RNA2}" \
                        --species "$SPECIES" --strain "$BASE" --cpus "$THREADS" --pasa_db sqlite
                else
                    funannotate train -i "$MASK" -o "$PREDICT" \
                        --single "${RNA_DIR}/${RNA1}" --species "$SPECIES" --strain "$BASE" \
                        --cpus "$THREADS"
                fi
                check_file "$PREDICT/training/funannotate_train.pasa.gff3" "Train output"
            else
                echo "Training output exists; skipping training."
            fi
        else
            echo "No RNA-seq data; skipping training."
        fi

        echo "--- Step 5: funannotate predict ---"
        [[ ! -s "$PREDICT_STATS" ]] && funannotate predict -i "$MASK" -o "$PREDICT" \
            --species "$SPECIES" --strain "$BASE" --cpus "$THREADS" --name "$BASE" || echo "Prediction exists; skipping."
        check_file "$PREDICT_STATS" "Prediction stats JSON"

        if (( RNAseq_FLAG )); then
            echo "--- Step 6:  funannotate update ---"
            unset PASACONF
            if [[ ! -s "$UPDATE" ]]; then
                funannotate update -i "$PREDICT" 
                check_file "$UPDATE" "Update output"
            else
                echo "Update output exists; skipping update."
            fi
        else
            echo "No RNA-seq data; skipping update."
        fi

        echo "--- Step 7a: eggnog mapper ---"
        [[ ! -s "$PREDICT/annotate_misc/eggnog.emapper.annotations" ]] && \
            MycoAssembler_eggnog.py -i "$PREDICT/predict_results/${SPECIES}_${BASE}.proteins.fa" \
                           -o "$PREDICT/annotate_misc/eggnog" --cpus "$THREADS" || \
            echo "Eggnog mapper output exists; skipping."
        check_file "$PREDICT/annotate_misc/eggnog.emapper.annotations" "Eggnog mapper annotations output"

        echo "--- Step 7b: funannotate annotate ---"
        [[ ! -s "$ANNO_OUT" ]] && funannotate annotate -i "$PREDICT" -o "$PREDICT" \
            --species "$SPECIES" --strain "$BASE" --cpus "$THREADS" || echo "Annotation exists; skipping."
        check_file "$ANNO_OUT" "Annotation stats JSON"

        echo "--- Completed annotation for $BASE ---"
    } > "$log_temp" 2>&1

    sed -E 's/\x1B\[[0-9;]*[mK]//g' "$log_temp" > "$log_file"
    rm -f "$log_temp"
}

run_busco() {
    local BASE="$1"        # sanitized DES
    local ASSEMBLER="$2"   # flye, canu, spades
    local SPECIES="$3"
    local BUSCO_LIN="$4"
    local MODE="$5"        # protein, genome, or both

    BASE="${BASE//[^a-zA-Z0-9._-]/_}"

    local PROT GEN
    PROT="$ANNO_DIR/$BASE.$ASSEMBLER/$BASE.funannotate04/annotate_results/${SPECIES}_${BASE}.proteins.fa"
    GEN="$ANNO_DIR/$BASE.$ASSEMBLER/$BASE.funannotate04/annotate_results/${SPECIES}_${BASE}.contigs.fsa"

    local timestamp log_temp log_file
    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part5.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part5.$BASE.$timestamp.log"

    {
        echo "--- BUSCO for $BASE ($ASSEMBLER) ---"
        echo "Threads: $THREADS"
        echo "Lineage: $BUSCO_LIN | Mode: $MODE"
        echo

        if [[ "$MODE" == "protein" || "$MODE" == "both" ]]; then
            if [[ -s "$PROT" ]]; then
                echo "Running BUSCO (protein) on $PROT"
                busco -i "$PROT" \
                      --out_path "$BUSCO_DIR/$BASE.$ASSEMBLER.$BUSCO_LIN.protein.buscoout" \
                      -l "$BUSCO_LIN" -m protein -c "$THREADS"
            else
                echo "Warning: Protein file not found: $PROT"
            fi
        fi

        if [[ "$MODE" == "genome" || "$MODE" == "both" ]]; then
            if [[ -s "$GEN" ]]; then
                echo "Running BUSCO (genome) on $GEN"
                busco -i "$GEN" \
                      --out_path "$BUSCO_DIR/$BASE.$ASSEMBLER.$BUSCO_LIN.genome.buscoout" \
                      -l "$BUSCO_LIN" -m genome -c "$THREADS"
            else
                echo "Warning: Genome file not found: $GEN"
            fi
        fi

        echo "--- BUSCO completed for $BASE ---"
    } > "$log_temp" 2>&1

    sed -E 's/\x1B\[[0-9;]*[mK]//g' "$log_temp" > "$log_file"
    rm -f "$log_temp"
}
