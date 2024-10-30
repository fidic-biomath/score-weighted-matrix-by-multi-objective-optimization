#!/bin/bash

for p1 in fit9.2 ; do #fit9.1; do
    for p4 in fit3.2 ; do #fit9.1;do
        for p9 in rmsd3.2 ; do #fit9.2 clustermsd3;do
            for p6 in rmsd3.2 ; do #fit9.2;do
                for p7 in rmsd3.2 ; do #fit9.1;do
                    for p2 in clustermsd3 rmsd9.1 fit9.1 rmsd3.1 fit9.2 rmsd9.2 clustermsd9; do
                        for p3 in clusterfit9 clustermsd9 fit9.1 fit3.1;do
                            for p5 in rmsd3.1 rmsd3.2 fit9.1;do
                                for p8 in fit9.2 clustermsd9 fit9.1 rmsd3.1 rmsd9.1;do
                                    echo "A C D E F G H I K L M N P Q R S T V W Y"> ${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==2 {print $0}'  ../matrices/${p1} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==3 {print $0}'  ../matrices/${p2} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==4 {print $0}'  ../matrices/${p3} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==5 {print $0}'  ../matrices/${p4} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==6 {print $0}'  ../matrices/${p5} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==7 {print $0}'  ../matrices/${p6} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==8 {print $0}'  ../matrices/${p7} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==9 {print $0}'  ../matrices/${p8} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                    awk 'NR==10 {print $0}' ../matrices/${p9} >>${p1}-${p2}-${p3}-${p4}-${p5}-${p6}-${p7}-${p8}-${p9}
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
