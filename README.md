// README
// =======================================
// Author: Olga Nikolova
// E-mail: olga.nikolova@gmail.com
// =======================================
// Preprocessing gene fusion data
// Input directory:
// Final output directory:
// 
// Preprocessing pipeline (ready to run)
//
// ==================================================
// Step 1: Check all unique file names present in data
// Output: list of file names in "./d2db_files_present.log" (n=4218)
// ==================================================
./check_data_d2db.sh

// ==================================================
// Step 2: Check the mapping file  from experiments to TCGA barcodes
// Output: 
// 1) List of files in data which are also mapped to TCGA: "db2d_files_present.log" (n=3972)
// 2) List of files in mapping file (experimental id mapped to TCGA) which however were not present in our data "db2d_files_missing.log" (n=3954)
// ==================================================
./check_data_d2db.sh

// This leaves n=246 experiments that were not mapped to TCGA

// ==================================================
// Process files: see file for comments
// Output: in mapped2TCGA_final
// ==================================================
./run.sh

// ==================================================
// Logs generated
// ==================================================
d2db_files_present.log
db2d_files_present.log
db2d_files_missing.log

// ==================================================
// Scripts
// ==================================================
./check_data_d2db.sh
./check_data_d2db.sh
./run.sh
