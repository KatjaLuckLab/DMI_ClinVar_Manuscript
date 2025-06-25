#! /bin/bash

python3 /AlphaFold/scripts/organize_computed_msas.py -computed_msas_path /AlphaFold/computed_msas/ -run_path /AlphaFold/HuRI_DMI_AF_modelling/AF_predictions/

NVIDIA_VISIBLE_DEVICES=1 time singularity run --contain --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \
--fasta_paths=/AlphaFold/HuRI_DMI_AF_modelling/AF_predictions/DMI.fasta \
--output_dir=/AlphaFold/HuRI_DMI_AF_modelling/AF_predictions/ \
--model_preset=multimer \
--db_preset=full_dbs \
--max_template_date=2020-05-14 \
--num_multimer_predictions_per_model=1 \
--use_gpu_relax=True \
--data_dir=/mnt/storage/alphafold/v232 \
--bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
--mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \
--obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \
--pdb_seqres_database_path=/mnt/storage/alphafold/v232/pdb_seqres/pdb_seqres.txt \
--template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \
--uniprot_database_path=/mnt/storage/alphafold/v232/uniprot/uniprot.fasta \
--uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \
--uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \
--use_precomputed_msas=True
