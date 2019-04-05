#!/bin/bash

python ATH_paser.py sg_all_gid.txt C > C_sg_all.txt
python ATH_paser.py sg_frameShift_gid.txt C > C_frameShift.txt
python ATH_paser.py sg_initiationAlt_gid.txt C > C_initiationAlt.txt
python ATH_paser.py sg_spliceSite_gid.txt C > C_spliceSite.txt
python ATH_paser.py sg_stopCodon_gid.txt C > C_stopCodon.txt
python ATH_paser.py sg_stopExtension_gid.txt C > C_stopExtension.txt

python ATH_paser.py sg_all_gid.txt F > F_sg_all.txt
python ATH_paser.py sg_frameShift_gid.txt F > F_frameShift.txt
python ATH_paser.py sg_initiationAlt_gid.txt F > F_initiationAlt.txt
python ATH_paser.py sg_spliceSite_gid.txt F > F_spliceSite.txt
python ATH_paser.py sg_stopCodon_gid.txt F > F_stopCodon.txt
python ATH_paser.py sg_stopExtension_gid.txt F > F_stopExtension.txt

python ATH_paser.py sg_all_gid.txt P > P_sg_all.txt
python ATH_paser.py sg_frameShift_gid.txt P > P_frameShift.txt
python ATH_paser.py sg_initiationAlt_gid.txt P > P_initiationAlt.txt
python ATH_paser.py sg_spliceSite_gid.txt P > P_spliceSite.txt
python ATH_paser.py sg_stopCodon_gid.txt P > P_stopCodon.txt
python ATH_paser.py sg_stopExtension_gid.txt P > P_stopExtension.txt

python fisher_test_modified.py C_sg_all.txt C_frameShift.txt > C_frameShift.out
python fisher_test_modified.py C_sg_all.txt C_initiationAlt.txt > C_initiationAlt.out
python fisher_test_modified.py C_sg_all.txt C_spliceSite.txt > C_spliceSite.out
python fisher_test_modified.py C_sg_all.txt C_stopCodon.txt > C_stopCodon.out
python fisher_test_modified.py C_sg_all.txt C_stopExtension.txt > C_stopExtension.out

python fisher_test_modified.py F_sg_all.txt F_frameShift.txt > F_frameShift.out
python fisher_test_modified.py F_sg_all.txt F_initiationAlt.txt > F_initiationAlt.out
python fisher_test_modified.py F_sg_all.txt F_spliceSite.txt > F_spliceSite.out
python fisher_test_modified.py F_sg_all.txt F_stopCodon.txt > F_stopCodon.out
python fisher_test_modified.py F_sg_all.txt F_stopExtension.txt > F_stopExtension.out

python fisher_test_modified.py P_sg_all.txt P_frameShift.txt > P_frameShift.out
python fisher_test_modified.py P_sg_all.txt P_initiationAlt.txt > P_initiationAlt.out
python fisher_test_modified.py P_sg_all.txt P_spliceSite.txt > P_spliceSite.out
python fisher_test_modified.py P_sg_all.txt P_stopCodon.txt > P_stopCodon.out
python fisher_test_modified.py P_sg_all.txt P_stopExtension.txt > P_stopExtension.out
