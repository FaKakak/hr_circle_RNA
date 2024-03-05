import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class Record_generator {
    public static final int LONG_CHAIN_LEN = 30;
    String file1 = "picture.txt";
    String file2 = "monitor.txt";
    Environment environment;

    public Record_generator(Environment environment) {
        this.environment = environment;
    }

    public void record(int h, int i) {
        int[] ch_num = new int[LONG_CHAIN_LEN], ch_nr_num = new int[LONG_CHAIN_LEN],
                ch_hr_num = new int[LONG_CHAIN_LEN], ch_0_num = new int[LONG_CHAIN_LEN],ch_hrnr_num = new int[LONG_CHAIN_LEN];
        int nr=0,cir_nr= 0,hr= 0,cir_hr= 0,hr_nr=0,cir_hr_nr= 0,ctr1= 0,ctr2= 0,ctr3= 0,
                cir_ctr1= 0,cir_ctr2= 0,cir_ctr3= 0,cir_unit= 0,total_mat_num= 0,unit= 0,raw_num= 0,
                long_chain_num, si;
        for (int y = 0; y < Environment.SIDE; y++) {
            for (int x = 0; x < Environment.SIDE; x++) {
                raw_num += environment.raw_arr[y][x];
                for (RNA p = environment.room_head[h][y][x].next;
                p != environment.room_head[h][y][x];
                p = p . next)
                {
                    unit++;
                    total_mat_num += p . length1 + p . length2;

                    int flag = findseq(RNA.NRSEQ, p);
                    int flag4 = findseq(RNA.HRSEQ, p);
                    if (flag == 0) nr++;
                    if (flag4 == 0) hr++;
                    if (flag == 0 && flag4 == 0) hr_nr++;

                    int flag1 = findseq(RNA.CTR1SEQ, p);
                    if (flag1 == 0) ctr1++; //
                    int flag2 = findseq(RNA.CTR2SEQ, p);
                    if (flag2 == 0) ctr2++; //
                    int flag3 = findseq(RNA.CTR3SEQ, p);
                    if (flag3 == 0) ctr3++; //

                    if (p . type1 == 1) {
                        cir_unit++;
                        if (flag == 0) cir_nr++;
                        if (flag4 == 0) cir_hr++;
                        if (flag == 0 && flag4 == 0) cir_hr_nr++;
                        if (flag1 == 0) cir_ctr1++;
                        if (flag2 == 0) cir_ctr2++;
                        if (flag3 == 0) cir_ctr3++;
                    }
                }
            }
        }
        total_mat_num += raw_num;

        System.out.printf("----- step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
                nr, cir_nr, hr, cir_hr, hr_nr, cir_hr_nr, ctr1, cir_ctr1, ctr2, cir_ctr2, ctr3, cir_ctr3,
                unit, cir_unit, total_mat_num, raw_num);
        try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
            writer.printf("----- step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
                    nr, cir_nr, hr, cir_hr, hr_nr, cir_hr_nr, ctr1, cir_ctr1, ctr2, cir_ctr2, ctr3, cir_ctr3,
                    unit, cir_unit, total_mat_num, raw_num);
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        }
        try (PrintWriter writer = new PrintWriter(new FileWriter(file2, true))) {
            writer.printf("----- step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
                    nr, cir_nr, hr, cir_hr, hr_nr, cir_hr_nr, ctr1, cir_ctr1, ctr2, cir_ctr2, ctr3, cir_ctr3,
                    unit, cir_unit, total_mat_num, raw_num);
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        }

        for (si = 0; si < LONG_CHAIN_LEN; si++)  //Ma2---start
        {
            ch_num[si] = 0;
            ch_nr_num[si] = 0;
            ch_hr_num[si] = 0;
            ch_hrnr_num[si] = 0;
            ch_0_num[si] = 0;
        }
        long_chain_num = 0;
        for (int y = 0; y < Environment.SIDE; y++) {
            for (int x = 0; x < Environment.SIDE; x++) {
                for (RNA p = environment.room_head[0][y][x].next;
                p != environment.room_head[0][y][x];
                p = p . next)
                {
                    ch_num[p . length1 - 1]++;
                    if (p . length1 > LONG_CHAIN_LEN) {
                        long_chain_num++;
                        for (int t = 0; t < p . length1; t++) {
                            char put=0;
                            switch (p.information[0][t]) {
                                case 'A' -> {
                                    put='A';
                                }
                                case 'C' -> {
                                    put='C';
                                }
                                case 'G' -> {
                                    put='G';
                                }
                                case 'U' -> {
                                    put = 'U';
                                }
                            }
                            System.out.printf(String.valueOf(put));
                            try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
                                writer.printf(String.valueOf(put));
                            } catch (IOException E) {
                                System.out.println("cannot open file");
                                System.exit(-1);
                            }
                        }
                        System.out.printf("\t%d\t%d\n", p . length1, p . type1);
                        try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
                            writer.printf("\t%d\t%d\n", p . length1, p . type1);
                        } catch (IOException E) {
                            System.out.println("cannot open file");
                            System.exit(-1);
                        }
                        for (int t = 0; t < p . length1; t++) {
                            char put=0;
                            switch (p.information[1][t]) {
                                case '0' -> {
                                    put = '-';
                                }
                                case 'A' -> {
                                    put = 'A';
                                }
                                case 'C' -> {
                                    put = 'C';
                                }
                                case 'G' -> {
                                    put = 'G';
                                }
                                case 'U' -> {
                                    put = 'U';
                                }
                            }
                            System.out.printf(String.valueOf(put));
                            try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
                                writer.printf(String.valueOf(put));
                            } catch (IOException E) {
                                System.out.println("cannot open file");
                                System.exit(-1);
                            }
                        }
                        System.out.printf("\t\t%d\n", p . type2);
                        try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
                            writer.printf("\t\t%d\n", p . type2);
                        } catch (IOException E) {
                            System.out.println("cannot open file");
                            System.exit(-1);
                        }
                    }
                    int flag1 = findseq(RNA.NRSEQ, p);
                    //flag3 = findseq(hrseq, hrlength, p);
                    int flag3 = findseq(RNA.HRSEQ, p);
                    if (flag1 == 0 && flag3 != 0) ch_nr_num[p . length1 - 1]++;
                    else if (flag1 != 0 && flag3 == 0) ch_hr_num[p . length1 - 1]++;
                    else if (flag1 == 0 && flag3 == 0) ch_hrnr_num[p . length1 - 1]++;
                    else ch_0_num[p . length1 - 1]++;
                }
            }
        }

        for (si = 0; si < LONG_CHAIN_LEN; si++) {
            System.out.printf("%dnt-%d(%d|%d/%d^%d), ", si + 1, ch_num[si], ch_hr_num[si], ch_nr_num[si], ch_hrnr_num[si], ch_0_num[si]);
            try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
                writer.printf("%dnt-%d(%d|%d/%d^%d), ", si + 1, ch_num[si], ch_hr_num[si], ch_nr_num[si], ch_hrnr_num[si], ch_0_num[si]);
            } catch (IOException E) {
                System.out.println("cannot open file");
                System.exit(-1);
            }
        }
        System.out.printf("\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);
        try (PrintWriter writer = new PrintWriter(new FileWriter(file1, true))) {
            writer.printf("\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        } //Ma2---end
    }


    private int findseq(char[] seq, RNA p) {
        int seqlength = seq.length;
        int length = RNA.MAX_RNA_LENGTH + RNA.MAX_CHAR_LENGTH;
        char[] inf = new char[length];    //RRRRR-ma3,  #define MAX_CHAR_LENGTH 20, --- to avoid using array out of bounds in the "dangerous block"
        //RRRRR-ma3,
        Arrays.fill(inf, '0');

        for (int a = 0; a < p.length1 + seqlength; a++)  // dangerous block
        {
            if (a < p.length1) inf[a] = p.information[0][a];
            else inf[a] = p.information[0][a - p.length1];
        }

        int flag2 = 0;
        // search for the subsequence
        if (p.length1 >= seqlength) {
            if (p.type1 == 0) {
                for (int b = 0; p.length1 - seqlength - b >= 0; b++) {
                    flag2 = 0;
                    for (int a = 0; a < seqlength; a++) {
                        if (inf[b + a] == seq[a]) continue;
                        else {
                            flag2 = 1;
                            break;
                        }  //this location has not this subsequence
                    }
                    if (flag2 == 0) break; // this location has this subsequence
                }
            } else if (p.type1 == 1) {

                for (int b = 0; b <= p.length1 - 1; b++) {
                    flag2 = 0;
                    for (int a = 0; a < seqlength; a++) {
                        if (inf[b + a] == seq[a]) continue;
                        else {
                            flag2 = 1;
                            break;
                        }
                    }
                    if (flag2 == 0) break;
                }
            }
        } else flag2 = 1;

        if (flag2 == 0) return (0);   //Yes, the sequence contains the subsequence
        else return (1);   //no, the sequence does not contain the subsequence
    }
}
