import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class Unit_case_generator {
    public static final double PSBP = 0.5;
    public static final double PBB = 0.00001;
    public static final double PLMC = 0.0000001;
    public static final double PEL = 0.0000001;
    public static final int MINI_CIRCLE_LENGTH = 8;
    public static final int FDA = 5;
    public static final int FHR = 1;
    public static final double FLT = 0.5;
    public static final double PAT = 0.5;
    public static final double PLT = 0.2;
    public static final double PMR = 0.01;
    public static final double PMV = 0.002;
    public static final double PMF = 0.001;
    public static final double PMFS = 0.9;
    public static final int TNSS = 1;
    public static final double PMD = 0.05;
    public static final double PNDE = 0.001;
    public static final double PFP = 0.01;

    public double RMRW(RNA p) {
        return pow(p.length1 + p.length2, 1 / 2.0);
    }

    Environment environment;
    int N;
    int[][] raw_arr;
    RNA[][][] cell_head;
    Random random = new Random(Main.SEED);
    int h;

    public Unit_case_generator(Environment environment) {
        this.environment = environment;
        N = Environment.SIDE;
        raw_arr = environment.raw_arr;
        cell_head = environment.room_head;
    }

    public void unit_case(int h) {
        this.h = h;

        List<Integer> xy_init = new LinkedList<>();
        for (int i = 0; i < Environment.CELLNUM; i++) {
            xy_init.add(i);
        }
        int xyIndex, xy, x, y;

        for (int i = 0; i < Environment.CELLNUM; i++) {
            // 随机选x,y，相当于xy_choose(void)
            xyIndex = random.nextInt(xy_init.size());
            xy = xy_init.get(xyIndex);
            xy_init.remove(xyIndex);
            x = xy % N;
            y = xy / N;

            raw(y, x);

            for (RNA p = cell_head[h][y][x].next; p != cell_head[h][y][x]; p = p.next) {
                switch (random.nextInt(6)) {
                    case 0 -> p = case0(p, y, x);
                    case 1 -> p = case1(p, y, x);
                    case 2 -> p = case2(p, y, x);
                    case 3 -> p = case3(p, y, x);
                    case 4 -> p = case4(p, y, x);
                    case 5 -> {
                        int[] rotate = case5(p, y, x);
                        p = fresh_unit(p, rotate[0], rotate[1]);
                    }
                }
            }
        }
    }

    private RNA fresh_unit(RNA p, int y, int x) {
        // 把当前RNA剥离这一层，搬到下一层，返回值为它的前一个RNA
        // 游离的RNA单独处理
        RNA p1 = p.prior;
        p.removeThis();
        cell_head[h ^ 1][y][x].addAfter(p);
        return p1;
    }

    private int findSeq(char[] seq, RNA p) {
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

    private void raw(int y, int x) {
        int raw_bef = raw_arr[y][x];     //Events of raw materials
        for (int k = 0; k < raw_bef; k++) {
            int randcaser = random.nextInt(2);
            switch (randcaser) {
                case 0 -> {  //forming nt
                    if (random.nextDouble() < PMF) {
                        raw_arr[y][x]--;
                        RNA p3 = new RNA();
                        int randnt = random.nextInt(4) + 1;
                        switch (randnt) {
                            case 1 -> p3.information[0][0] = 'A';
                            case 2 -> p3.information[0][0] = 'C';
                            case 3 -> p3.information[0][0] = 'G';
                            case 4 -> p3.information[0][0] = 'U';
                        }

                        p3.length1 = 1;
                        fresh_unit(p3, y, x);
                    }
                }
                case 1 -> {   // raw moving
                    if (random.nextDouble() < PMR)   //Ma2: with toroidal topology to avoid edge effects
                    {
                        raw_arr[y][x]--;   //Ma2
                        int randcaser1 = random.nextInt(4);   // Four possible directions
                        switch (randcaser1) {
                            case 0 -> raw_arr[y][(N + x - 1) % N]++;  //Na2: toroidal topology
                            case 1 -> raw_arr[y][(x + 1) % N]++;
                            case 2 -> raw_arr[(N + y - 1) % N][x]++;
                            case 3 -> raw_arr[(y + 1) % N][x]++;
                        }
                    }
                }
            }
        }
    }

    private RNA case0(RNA p, int y, int x) { // Chain ligation with mineral catalysis
        if (p.type1 == 0)
        {
            if (p.length1 >= MINI_CIRCLE_LENGTH)      //Ma2: considering self-ligation before ligation between RNAs //new
            {
                double f = PEL;
                int flag = 0;  //Ma2
                int hrlength = RNA.HRSEQ.length;
                for (int a = 0; a < hrlength / 2; a++)  //Ma2
                {
                    if (p.information[0][p.length1 - hrlength / 2 + a] == RNA.HRSEQ[a] && p.information[1][p.length1 - hrlength / 2 + a] == '0')continue;
                    else { flag = 1; break; }
                }
                if (flag == 0)
                {
                    for (int a = 0; a < hrlength / 2; a++)
                    {
                        if (p.information[0][a] == RNA.HRSEQ[hrlength / 2 + a] && p.information[1][a] == '0')continue;
                        else { flag = 1; break; }
                    }
                }
                if (flag == 0)f = f * FHR;  //Ma2
                if (random.nextDouble() < f)     //Ma2
                {
                    p.type1 = 1;
                    return fresh_unit(p,y,x);
                }
            }

            for (RNA p3 = p.next; p3 != p; p3 = p3.next)
            {
                if (p3 == environment.room_head[h][y][x]) { p3 = environment.room_head[h][y][x].next; if (p3 == p)break; }
                if (p3.type1 == 0)   //Ma2---start: the p3 needs not to be single chain
                {
                    if (random.nextDouble() < PLMC / (p.length1 * p3.length1))   //Ma2: the random ligation should be influenced by the lengths of both strands
                    {
                        if (p.length1 + p3.length1 > RNA.MAX_RNA_LENGTH - 1)
                        {
                            continue;
                        }

                        for (int a = 0; a < p3.length1; a++)
                        {
                            p.information[0][a + p.length1] = p3.information[0][a];
                            p.information[1][a + p.length1] = p3.information[1][a];  //Ma2
                        }

                        for (c2_frag c2f3 = p3.chain2.next; c2f3 != p3.chain2; c2f3 = c2f3.next)  //Ma2---start
                        {
                            c2f3.start += p.length1;
                        }
                        p.chain2.prior.next = p3.chain2.next;
                        p3.chain2.next.prior = p.chain2.prior;
                        p.chain2.prior = p3.chain2.prior;
                        p3.chain2.prior.next = p.chain2;   //Ma2---end

                        p.length1 = p.length1 + p3.length1;
                        p.length2 = p.length2 + p3.length2;  //Ma2

                        p3.removeThis();
                        break;
                    }
                }
            }
        }
        return fresh_unit(p,y,x);
    }

    private RNA case1(RNA p, int y, int x) {
        if (p.length1 == 1)  // Decay of mononucleotide
        {
            //if (p.information[0][0] == 0) { printf("unexpected error on p-nucleotide"); exit(0); }
            if (p.length2 == 0)
            {
                if (random.nextDouble() < PMD)
                {
                    raw_arr[y][x]++;
                    RNA rP = p.prior;
                    p.removeThis();
                    return rP;
                }
            }
            else if (p.length2 == 1)      //Ma2: The decay of the paired nucleotides at the same time
            {
                if (random.nextDouble() < PMD * sqrt(PMD))
                {
                    raw_arr[y][x] += 2;
                    RNA rP = p.prior;
                    p.removeThis();
                    return rP;
                }
            }
        }
        else                  //Degradation of chain
        {
            if (p.type1 == 0)
            {
                c2_frag c2f1 = p.chain2.prior;

                // Nucleotide residue decaying at the end of RNA
                if (p.information[1][p.length1 - 1] == '0')         // Single chain at the end
                {
                    if (random.nextDouble() < PNDE)
                    {
                        p.information[0][p.length1 - 1] = '0';
                        p.length1--;
                        raw_arr[y][x]++;
                    }
                }
                else if (p.information[1][p.length1 - 1] != '0')  //Double chain at the end
                {
                    if (random.nextDouble() < PNDE * sqrt(PNDE))                 // The decay of the paired residues at the same time
                    {
                        p.information[0][p.length1 - 1] = '0';
                        p.information[1][p.length1 - 1] = '0';
                        raw_arr[y][x]+=2;
                        p.length1--;
                        p.length2--;
                        c2f1.length--;

                        if (c2f1.length == 0)
                        {
                            c2f1.removeThis();
                        }
                    }
                }

                if (p.length1 == 1)
                {
                    return fresh_unit(p,y,x);
                }

                // Nucleotide residue decaying at the start of RNA
                c2_frag c2f2 = p.chain2.next;
                if (p.information[1][0] == '0')      //Single chain at the start
                {
                    if (random.nextDouble() < PNDE)             // The decay of the start residue on this single chain
                    {
                        for (int b = 1; b < p.length1; b++)
                        {
                            p.information[0][b - 1] = p.information[0][b];
                            p.information[1][b - 1] = p.information[1][b];
                        }
                        p.information[0][p.length1 - 1] = '0';
                        p.information[1][p.length1 - 1] = '0';
                        p.length1--;
                        raw_arr[y][x]++;

                        for (c2_frag c2f3 = c2f2; c2f3 != p.chain2; c2f3 = c2f3.next)    //Ma2
                        {
                            c2f3.start--;
                        }
                    }
                }
                else if (p.information[1][0] != '0') //Double chain at the start
                {
                    if (random.nextDouble() < PNDE * sqrt(PNDE)) { // The decay of the paired residues at the same time

                        for (int b = 1; b < p.length1; b++)
                        {
                            p.information[0][b - 1] = p.information[0][b];
                            p.information[1][b - 1] = p.information[1][b];
                        }
                        p.information[0][p.length1 - 1] = '0';
                        p.information[1][p.length1 - 1] = '0';
                        raw_arr[y][x]+=2;

                        p.length1--;
                        p.length2--;
                        c2f2.length--;

                        for (c2_frag c2f3 = c2f2.next; c2f3 != p.chain2; c2f3 = c2f3.next)
                        {
                            c2f3.start--;
                        }

                        if (c2f2.length == 0)
                        {
                            c2f2.removeThis();
                        }
                    }
                }

                if (p.length1 == 1)
                {
                    return fresh_unit(p,y,x);
                }

                while (true)
                {
                    int j;
                    int hrlength = RNA.HRSEQ.length;
                    for (j = p.length1; j > 1; j--)
                    {
                        double f = PBB;
                        for (c2f1 = p.chain2.prior; c2f1 != p.chain2; c2f1 = c2f1.prior)
                        {
                            if (j > c2f1.start + 1)
                            {
                                if (j <= c2f1.length + c2f1.start)         // Falling into double chain region
                                {
                                    f = f * sqrt(f);
                                }
                                break;
                            }
                        }

                        int flag = 1;
                        if (j - 1 >= hrlength / 2 && p.length1 - j + 1 >= hrlength / 2)       //new
                        {
                            flag = 0;
                            for (int a = 0; a < hrlength; a++)
                            {
                                if (p.information[0][j - 1 - hrlength / 2 + a] == RNA.HRSEQ[a] && p.information[1][j - 1 - hrlength / 2 + a] == '0')continue;
                                else { flag = 1; break; }
                            }
                        }
                        if (flag == 0)f = f * FHR;   //Ma2

                        c2f1 = p.chain2.prior;
                        c2f2 = p.chain2.next;
                        if (random.nextDouble() < f)
                        {
                            RNA p3 = new RNA();
                            for (int b = 0; b < p.length1 - j + 1; b++)
                            {
                                p3.information[0][b] = p.information[0][b + j - 1];
                                p.information[0][b + j - 1] = '0';
                            }
                            p3.length1 = p.length1 - j + 1;
                            p.length1 = j - 1;
                            if (c2f1 != p.chain2 && c2f1.start + c2f1.length > j - 1)   //Ma2
                            {
                                for (int b = 0; b < c2f1.start + c2f1.length - j + 1; b++)
                                {
                                    p3.information[1][b] = p.information[1][b + j - 1];
                                    p.information[1][b + j - 1] = '0';
                                }

                                while (c2f1 != p.chain2)
                                {
                                    c2_frag c2f3 = new c2_frag();
                                    c2f3.prior = p3.chain2;
                                    c2f3.next = p3.chain2.next;
                                    p3.chain2.next.prior = c2f3;
                                    p3.chain2.next = c2f3;

                                    if (j <= c2f1.start + 1)
                                    {
                                        c2f3.length = c2f1.length;
                                        c2f3.start = c2f1.start - j + 1;
                                        p3.length2 = p3.length2 + c2f3.length;
                                        p.length2 = p.length2 - c2f3.length;

                                        (c2f1.prior).next = c2f1.next;
                                        (c2f1.next).prior = c2f1.prior;
                                        c2f1 = c2f1.prior;

                                        if (j >= c2f1.start + c2f1.length + 1)
                                        {
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        c2f3.length = c2f1.start + c2f1.length - j + 1;
                                        c2f3.start = 0;
                                        c2f1.length = c2f1.length - c2f3.length;

                                        p3.length2 = p3.length2 + c2f3.length;
                                        p.length2 = p.length2 - c2f3.length;
                                        break;
                                    }
                                }
                            }
                            fresh_unit(p3,y,x);
                            break;     //Bond break occurs
                        }
                    }
                    if (j == 1) break;
                }
            }
            else if (p.type1 == 1)
            {
                int flag5 = 0;  //Ma2: a flag labeling whether the circle chain breaking has happened
                for (int j = p.length1; j > 0; j--)
                {
                    double f = PBB;
                    int flag1 = 0; //Ma2
                    int hrlength = RNA.HRSEQ.length;
                    c2_frag c2f1;
                    for (c2f1 = p.chain2.prior; c2f1 != p.chain2; c2f1 = c2f1.prior)
                    {
                        if (j > c2f1.start + 1)
                        {
                            if (j <= c2f1.length + c2f1.start)         // Falling into double chain region
                            {
                                f = f * sqrt(f); flag1 = 1;  //Ma2
                            }
                            break;
                        }
                    }
                    if (j == 1 && p.type2 == 1)f = f * sqrt(f);

                    int flag,a;
                    for (flag = 0, a = 0; a < hrlength; a++)                    //new
                    {
                        if (j - 1 - hrlength / 2 + a >= p.length1)
                        {
                            if (p.information[0][j - 1 - hrlength / 2 + a - p.length1] == RNA.HRSEQ[a] && p.information[1][j - 1 - hrlength / 2 + a - p.length1] == '0')continue;
                            else { flag = 1; break; }
                        }
                        else if (j - 1 - hrlength / 2 + a < 0)
                        {
                            if (p.information[0][p.length1 + j - 1 - hrlength / 2 + a] == RNA.HRSEQ[a] && p.information[1][p.length1 + j - 1 - hrlength / 2 + a] == '0')continue;
                            else { flag = 1; break; }
                        }
                        else
                        {
                            if (p.information[0][j - 1 - hrlength / 2 + a] == RNA.HRSEQ[a] && p.information[1][j - 1 - hrlength / 2 + a] == '0')continue;
                            else { flag = 1; break; }
                        }
                    }
                    if (flag == 0)f = f * FHR;   //Ma2

                    if (random.nextDouble() < f)
                    {
                        flag5 = j;  // Ma2: the location of the circle chain breaking
                        p.type1 = 0;
                        p.type2 = 0;
                        if (j == 1)break;

                        char[][] temp_information = new char[2][RNA.MAX_RNA_LENGTH];
                        for (int b = 0; b < j - 1; b++)
                        {
                            temp_information[0][b] = p.information[0][b];
                            temp_information[1][b] = p.information[1][b];
                        }
                        for (int b = 0; b < p.length1; b++)
                        {
                            if (b < p.length1 - j + 1)
                            {
                                p.information[0][b] = p.information[0][b + j - 1];
                                p.information[1][b] = p.information[1][b + j - 1];
                            }
                            else
                            {
                                p.information[0][b] = temp_information[0][b - p.length1 + j - 1];
                                p.information[1][b] = temp_information[1][b - p.length1 + j - 1];
                            }
                        }
                        if (p.chain2.next == p.chain2)break;  //Ma2

                        c2_frag c2f3;
                        if (flag1 == 1) //f == PBB * sqrt(PBB)) //Ma2
                        {
                            c2f3 = new c2_frag();
                            c2f3.prior = c2f1;
                            c2f3.next = c2f1.next;
                            c2f1.next.prior = c2f3;
                            c2f1.next = c2f3;

                            c2f3.start = j - 1;
                            c2f3.length = c2f1.start + c2f1.length - j + 1;
                            c2f1.length = c2f1.length - c2f3.length;
                        }
                        else
                        {
                            c2f3 = c2f1.next;   //Ma2
                        }

                        //Ma2: two "if" sentences were deleted

                        for (c2_frag c2f2 = c2f3; c2f2 != p.chain2; c2f2 = c2f2.next) c2f2.start = c2f2.start - j + 1;
                        for (c2_frag c2f2 = p.chain2.next; c2f2 != c2f1.next; c2f2 = c2f2.next) c2f2.start = c2f2.start + p.length1 - j + 1;

                        //if (p.chain2.next.next == p.chain2)break;  //Ma2

                        if (c2f3 != p.chain2 && c2f1 != p.chain2)   //Ma2
                        {
                            p.chain2.prior.next = p.chain2.next;
                            p.chain2.next.prior = p.chain2.prior;
                            p.chain2.next = c2f3;
                            p.chain2.prior = c2f1;
                            c2f1.next = p.chain2;
                            c2f3.prior = p.chain2;
                        }

                        break;
                    }
                }
                if (flag5 != 0) //Ma2--start: Circlar chain breaking has happened
                {
                    int flag4 = p.length1 - flag5 + 1;
                    int hrlength = RNA.HRSEQ.length;
                    while (true)
                    {
                        int j;
                        for (j = p.length1; j > flag4; j--)  //Ma2: Only consider those bonds that have not been considered before the circular chain breaking
                        {
                            double f = PBB;
                            for (c2_frag c2f1 = p.chain2.prior; c2f1 != p.chain2; c2f1 = c2f1.prior)
                            {
                                if (j > c2f1.start + 1)
                                {
                                    if (j <= c2f1.length + c2f1.start)         // Falling into double chain region
                                    {
                                        f = f * sqrt(f);
                                    }
                                    break;
                                }
                            }

                            int flag = 1;
                            if (j - 1 >= hrlength / 2 && p.length1 - j + 1 >= hrlength / 2)       //new
                            {
                                flag = 0;
                                for (int a = 0; a < hrlength; a++)
                                {
                                    if (p.information[0][j - 1 - hrlength / 2 + a] == RNA.HRSEQ[a] && p.information[1][j - 1 - hrlength / 2 + a] == '0')continue;
                                    else { flag = 1; break; }
                                }
                            }
                            if (flag == 0)f = f * FHR;   //Ma2

                            c2_frag c2f1 = p.chain2.prior;
                            c2_frag c2f2 = p.chain2.next;
                            if (random.nextDouble() < f)
                            {
                                RNA p3 = new RNA();

                                for (int b = 0; b < p.length1 - j + 1; b++)
                                {
                                    p3.information[0][b] = p.information[0][b + j - 1];
                                    p.information[0][b + j - 1] = '0';
                                }
                                p3.length1 = p.length1 - j + 1;
                                p.length1 = j - 1;
                                if (c2f1 != p.chain2 && c2f1.start + c2f1.length > j - 1)   //Ma2
                                {
                                    for (int b = 0; b < c2f1.start + c2f1.length - j + 1; b++)
                                    {
                                        p3.information[1][b] = p.information[1][b + j - 1];
                                        p.information[1][b + j - 1] = '0';
                                    }

                                    while (c2f1 != p.chain2)
                                    {
                                        c2_frag c2f3 = new c2_frag();
                                        c2f3.prior = p3.chain2;
                                        c2f3.next = p3.chain2.next;
                                        p3.chain2.next.prior = c2f3;
                                        p3.chain2.next = c2f3;

                                        if (j <= c2f1.start + 1)
                                        {
                                            c2f3.length = c2f1.length;
                                            c2f3.start = c2f1.start - j + 1;
                                            p3.length2 = p3.length2 + c2f3.length;
                                            p.length2 = p.length2 - c2f3.length;

                                            (c2f1.prior).next = c2f1.next;
                                            (c2f1.next).prior = c2f1.prior;
                                            c2f1 = c2f1.prior;
                                            if (j >= c2f1.start + c2f1.length + 1)
                                            {
                                                break;
                                            }
                                        }
                                        else
                                        {
                                            c2f3.length = c2f1.start + c2f1.length - j + 1;
                                            c2f3.start = 0;
                                            c2f1.length = c2f1.length - c2f3.length;
                                            p3.length2 = p3.length2 + c2f3.length;
                                            p.length2 = p.length2 - c2f3.length;
                                            break;
                                        }
                                    }
                                }
                                fresh_unit(p3,y,x);
                                break;    //Bond break occurs
                            }
                        }
                        if (j == flag4) break;
                    }
                }//Ma2--end
            }
        }
        return fresh_unit(p,y,x);
    }

    private RNA case2(RNA p, int y, int x) {
        for (c2_frag c2f2 = p.chain2.next; c2f2 != p.chain2; c2f2 = c2f2.next)         //Template-directed ligation
        {
            double rtdaddlig = random.nextDouble();
            c2_frag c2f3 = c2f2.next;
            if (c2f3 != p.chain2)
            {
                if (c2f2.length + c2f2.start == c2f3.start && rtdaddlig < PLT)
                {
                    c2f2.length = c2f2.length + c2f3.length;
                    c2f3.removeThis();
                    break;
                }
            }
            else
            {
                if (p.type1 == 1 && c2f2.length + c2f2.start == p.length1 && p.chain2.next.start == 0 && rtdaddlig < PLT)
                {
                    //if (p.information[1][0] == 0 && p.information[1][p.length1 - 1] == 0) { printf("unexpected error on add-nucleotide2"); exit(0); }
                    c2_frag c2f1 = p.chain2.prior;
                    c2_frag c2f4 = p.chain2.next;   //Ma2
                    if (c2f1 == c2f4)
                    {
                        p.type2 = 1;
                        break;
                    }

                    char[][] temp_information = new char[2][RNA.MAX_RNA_LENGTH];
                    for (int b = 0; b < c2f1.length; b++)
                    {
                        temp_information[0][b] = p.information[0][b + p.length1 - c2f1.length];
                        temp_information[1][b] = p.information[1][b + p.length1 - c2f1.length];
                    }
                    for (int b = p.length1 - 1; b >= 0; b--)
                    {
                        if (b >= c2f1.length)
                        {
                            p.information[0][b] = p.information[0][b - c2f1.length];
                            p.information[1][b] = p.information[1][b - c2f1.length];
                        }
                        else
                        {
                            p.information[0][b] = temp_information[0][b];
                            p.information[1][b] = temp_information[1][b];
                        }
                    }

                    c2f4.length = c2f4.length + c2f1.length;   //Ma2
                    for (c2f4 = c2f4.next; c2f4 != p.chain2; c2f4 = c2f4.next)  //Ma2
                    {
                        c2f4.start = c2f4.start + c2f1.length;
                    }
                    c2f1.prior.next = p.chain2;
                    p.chain2.prior = c2f1.prior;
                    break;
                }
            }
        }

//if (p.type1 == 0) { fresh_unit(); break; }  //Ma2-1  Linear RNA cannot attract substrates
//Template-directed attraction of substrates
        double f1;
        if (p.type1 == 1) f1=PAT;  //Ma2-1
        else f1 = PAT * FLT;   //Ma2-1
        for (RNA p3 = p.next; p3 != p; p3 = p3.next)
        {
            if (p3 == environment.room_head[h][y][x])
            {
                p3 = environment.room_head[h][y][x].next;
                if (p3 == p)break;
            }
            if (p3.length2 == 0 && p3.length1 <= p.length1)
            {
                c2_frag c2f2 = p.chain2.next;
                if (p.type1 == 0 && (p.length2 == 0 || c2f2.start != 0))
                {
                    int length=0;
                    if (p.length2 == 0)
                    {
                        length = p.length1;
                    }
                    else if (c2f2.start != 0)
                    {
                        length = c2f2.start;
                    }

                    if (p3.length1 <= length)
                    {
                        int flag=1,c;
                        for (c = 0; c <= length - p3.length1; c++)
                        {
                            int b;
                            for (flag = 0, b = 0; b < p3.length1; b++)
                            {
                                if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b]) == 'G'+'C'))continue;
                                else if (random.nextDouble() < PFP)continue;
                                else { flag = 1; break; }
                            }
                            if (flag == 0)break;
                        }
                        if (flag == 0)
                        {
                            double rtdaddphili = random.nextDouble() * FDA;  //Ma2
                            if (rtdaddphili < f1)
                            {
                                for (int a = 0; a < p3.length1; a++)
                                {
                                    p.information[1][c + a] = p3.information[0][p3.length1 - 1 - a];
                                }

                                c2_frag c2f3 = new c2_frag();
                                c2f3.prior = p.chain2;
                                c2f3.next = p.chain2.next;
                                p.chain2.next.prior = c2f3;
                                p.chain2.next = c2f3;

                                c2f3.start = c;
                                c2f3.length = p3.length1;

                                p.length2 = p.length2 + p3.length1;
                                p3.removeThis();
                                break;
                            }
                        }
                    }
                }
                else if (p.type1 == 1 && p.length2 == 0)
                {
                    if (p.length1 >= p3.length1)
                    {
                        int flag=1;
                        int c;
                        for (c = 0; c <= p.length1 - 1; c++)
                        {
                            int b;
                            for (flag = 0, b = 0; b < p3.length1; b++)
                            {
                                if (c + b < p.length1)
                                {
                                    if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b]) == 'G'+'C'))continue;
                                    else if (random.nextDouble() < PFP)continue;
                                    else { flag = 1; break; }
                                }
                                else
                                {
                                    if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b - p.length1]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + b - p.length1]) == 'C'+'G'))continue;
                                    else if (random.nextDouble() < PFP)continue;
                                    else { flag = 1; break; }
                                }
                            }
                            if (flag == 0)break;
                        }
                        if (flag == 0)
                        {
                            char[][] temp_information = new char[2][RNA.MAX_RNA_LENGTH];
                            double rtdaddphili = random.nextDouble() * FDA;   //Ma2
                            if (rtdaddphili < f1)
                            {
                                c2_frag c2f3 = new c2_frag();
                                c2f3.prior = p.chain2;
                                c2f3.next = p.chain2.next;
                                p.chain2.next.prior = c2f3;
                                p.chain2.next = c2f3;
                                c2f3.length = p3.length1;

                                if (p3.length1 + c <= p.length1)
                                {
                                    for (int a = 0; a < p3.length1; a++) {
                                        p.information[1][c + a] = p3.information[0][p3.length1 - 1 - a];
                                    }

                                    c2f3.start = c;
                                }
                                else
                                {
                                    if (p.length1 - c >= 0)
                                        System.arraycopy(p.information[0], c, temp_information[0], 0, p.length1 - c);

                                    for (int b = p.length1 - 1; b >= 0; b--)
                                    {
                                        if (b >= p.length1 - c) {
                                            p.information[0][b] = p.information[0][b - p.length1 + c];
                                        }
                                        else {
                                            p.information[0][b] = temp_information[0][b];
                                        }
                                    }
                                    for (int a = 0; a < p3.length1; a++) { p.information[1][a] = p3.information[0][p3.length1 - 1 - a]; }

                                    c2f3.start = 0;
                                }
                                p.length2 = p.length2 + p3.length1;
                                p3.removeThis();
                                break;
                            }
                        }
                    }
                }
                else if (p.type1 == 1 && c2f2.start != 0)
                {
                    c2_frag c2f1 = p.chain2.prior;
                    int length = c2f2.start + p.length1 - c2f1.start - c2f1.length;

                    if (length >= p3.length1)
                    {
                        int flag = 1;
                        int c;
                        for (c = 0; c <= length - p3.length1; c++)
                        {
                            int b;
                            for (flag = 0, b = 0; b < p3.length1; b++)
                            {
                                if (c + c2f1.start + c2f1.length + b < p.length1)   //Ma2
                                {
                                    if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + c2f1.start + c2f1.length + b]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + c2f1.start + c2f1.length + b]) == 'G'+'C'))continue;  //Ma2
                                    else if (random.nextDouble() < PFP)continue;
                                    else { flag = 1; break; }
                                }
                                else
                                {
                                    if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + c2f1.start + c2f1.length + b - p.length1]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c + c2f1.start + c2f1.length + b - p.length1]) == 'G'+'C'))continue;  //Ma2
                                    else if (random.nextDouble() < PFP)continue;
                                    else { flag = 1; break; }
                                }
                            }
                            if (flag == 0)break;
                        }
                        if (flag == 0)
                        {
                            double rtdaddphili = random.nextDouble();
                            if (c != 0)rtdaddphili = rtdaddphili * FDA;  //Ma2
                            if (rtdaddphili < f1)
                            {
                                c2_frag c2f3 = new c2_frag();
                                c2f3.length = p3.length1;

                                if (p3.length1 + c <= p.length1 - c2f1.start - c2f1.length)
                                {
                                    for (int a = 0; a < p3.length1; a++)
                                        p.information[1][c2f1.start + c2f1.length + c + a] = p3.information[0][p3.length1 - 1 - a];

                                    c2f3.next = p.chain2;
                                    c2f3.prior = p.chain2.prior;
                                    p.chain2.prior.next = c2f3;
                                    p.chain2.prior = c2f3;
                                    c2f3.start = c + c2f1.start + c2f1.length;
                                }
                                else if (c >= p.length1 - c2f1.start - c2f1.length)
                                {
                                    for (int a = 0; a < p3.length1; a++)
                                        p.information[1][c - p.length1 + c2f1.start + c2f1.length + a] = p3.information[0][p3.length1 - 1 - a];

                                    c2f3.prior = p.chain2;
                                    c2f3.next = p.chain2.next;
                                    p.chain2.next.prior = c2f3;
                                    p.chain2.next = c2f3;
                                    c2f3.start = c - p.length1 + c2f1.start + c2f1.length;
                                }
                                else
                                {
                                    char[][] temp_information = new char[2][RNA.MAX_RNA_LENGTH];
                                    for (int b = 0; b < p.length1 - c2f1.start - c2f1.length - c; b++)
                                    {
                                        temp_information[0][b] = p.information[0][b + c2f1.start + c2f1.length + c];
                                        temp_information[1][b] = p.information[1][b + c2f1.start + c2f1.length + c];
                                    }
                                    for (int b = p.length1 - 1; b >= 0; b--)
                                    {
                                        if (b >= p.length1 - c2f1.start - c2f1.length - c)
                                        {
                                            p.information[0][b] = p.information[0][b - p.length1 + c2f1.start + c2f1.length + c];
                                            p.information[1][b] = p.information[1][b - p.length1 + c2f1.start + c2f1.length + c];
                                        }
                                        else
                                        {
                                            p.information[0][b] = temp_information[0][b];
                                            p.information[1][b] = temp_information[1][b];
                                        }
                                    }
                                    for (int a = 0; a < p3.length1; a++) { p.information[1][a] = p3.information[0][p3.length1 - 1 - a]; }

                                    for (c2f2 = p.chain2.next; c2f2 != p.chain2; c2f2 = c2f2.next)
                                        c2f2.start = c2f2.start + p.length1 - c2f1.start - c2f1.length - c;

                                    c2f3.prior = p.chain2;
                                    c2f3.next = p.chain2.next;
                                    p.chain2.next.prior = c2f3;
                                    p.chain2.next = c2f3;
                                    c2f3.start = 0;
                                }
                                p.length2 = p.length2 + p3.length1;

                                p3.removeThis();
                                break;
                            }
                        }
                    }
                }
                int flag3 = 0;  //Ma2:  a flag labeling whether the attraction has happened
                for (c2f2 = p.chain2.next; c2f2 != p.chain2; c2f2 = c2f2.next)
                {
                    int length;
                    if (c2f2.next == p.chain2)
                    {
                        if (p.type1 == 1 && p.chain2.next.start != 0) break;  //Ma2:  "p.length2 == 0 ||" is deleted
                        length = p.length1 - c2f2.start - c2f2.length;
                    }
                    else
                    {
                        length = c2f2.next.start - c2f2.start - c2f2.length;
                    }

                    if (p3.length1 <= length)
                    {
                        int flag = 1;
                        int c;
                        for (c = 0; c <= length - p3.length1; c++)
                        {
                            int b;
                            for (flag = 0, b = 0; b < p3.length1; b++)
                            {
                                if (((p3.information[0][p3.length1 - 1 - b] + p.information[0][c2f2.start + c2f2.length + c + b]) == 'A'+'U')||((p3.information[0][p3.length1 - 1 - b] + p.information[0][c2f2.start + c2f2.length + c + b]) == 'G'+'C'))continue;
                                else if (random.nextDouble() < PFP)continue;
                                else { flag = 1; break; }
                            }
                            if (flag == 0)break;
                        }
                        if (flag == 0)
                        {
                            double rtdaddphili = random.nextDouble();
                            if (c != 0)rtdaddphili = rtdaddphili * FDA; //Ma2
                            if (rtdaddphili < f1)
                            {
                                for (int a = 0; a < p3.length1; a++) {
                                    p.information[1][c2f2.start + c2f2.length + c + a] = p3.information[0][p3.length1 - 1 - a];
                                }

                                c2_frag c2f3 = new c2_frag();
                                c2f3.prior = c2f2;
                                c2f3.next = c2f2.next;
                                c2f2.next.prior = c2f3;
                                c2f2.next = c2f3;

                                c2f3.start = c2f2.start + c2f2.length + c;
                                c2f3.length = p3.length1;

                                p.length2 = p.length2 + p3.length1;

                                p3.removeThis();
                                flag3 = 1;   //Ma2:  the attraction has happened.  -- "flag=2" is changed.
                                break;
                            }
                        }
                    }
                }
                if (flag3 == 1) break;  //Ma2 "flag==2" is changed
            }
        }
        return fresh_unit(p,y,x);
    }

    private RNA case3(RNA p, int y, int x) { // Separation
        if (p.length2 != 0)
        {
            int j = 0;
            c2_frag c2f2;
            for (c2f2 = p.chain2.next; c2f2 != p.chain2; c2f2 = c2f2.next) j++;
            int randseq = random.nextInt(j);
            c2f2 = p.chain2.next;
            for (j = 1; j <= randseq; j++)
            {
                c2f2 = c2f2.next;
            }

            if (random.nextDouble() < pow(PSBP, sqrt(c2f2.length *1.0)))  //Ma2-1
            {
                RNA p3 = new RNA();
                for (int b = 0; b < c2f2.length; b++)
                {
                    p3.information[0][b] = p.information[1][c2f2.start + c2f2.length - 1 - b];
                    p.information[1][c2f2.start + c2f2.length - 1 - b] = '0';
                }
                p3.length1 = c2f2.length;
                p.length2 = p.length2 - c2f2.length;
                p3.length2 = 0;
                if (p.type2 == 1)
                {
                    p3.type1 = 1;
                    p.type2 = 0;    //Ma2
                }
                else p3.type1 = 0;
                p3.type2 = 0;      //Ma2
                c2f2.removeThis();
                fresh_unit(p3,y,x);
            }
        }
        return fresh_unit(p,y,x);
    }

    private RNA case4(RNA p, int y, int x) { // nr catalyses the synthesis of nt
        if (p.type1 == 0 && p.length2 == 0 && p.length1 >= RNA.NRSEQ.length && p.length1 < 2 * RNA.NRSEQ.length)  //Ma2-1, p.length1 < 2 * nrlength:: Nr cannot be much longer than its characteristic domain
        {
            int flag = findSeq(RNA.NRSEQ,p);
            if (flag == 0)     //Ma2-1
            {
                int nt_turn = TNSS;
                int raw_bef = raw_arr[y][x];
                for (int k = 0; k < raw_bef; k++)
                {
                    if (nt_turn <= 0)break;
                    nt_turn--;
                    if (random.nextDouble() < PMFS)
                    {
                        raw_arr[y][x]--;
                        RNA p3 = new RNA();
                        int randnt = random.nextInt(4) + 1;
                        switch (randnt) {
                            case 1 -> p3.information[0][0] = 'A';
                            case 2 -> p3.information[0][0] = 'C';
                            case 3 -> p3.information[0][0] = 'G';
                            case 4 -> p3.information[0][0] = 'U';
                        }

                        p3.length1 = 1;
                        fresh_unit(p3,y,x);
                    }
                }
            }
        }
        return fresh_unit(p,y,x);
    }

    private int[] case5(RNA p, int y, int x) {
        int ny = y, nx = x;
        if (random.nextDouble() * RMRW(p) < PMV)    //Ma2: with toroidal topology to avoid edge effects
        {
            int randcase1 = random.nextInt(4);   // Four possible directions
            switch (randcase1) {
                case 0 -> nx = (N + x - 1) % N;
                case 1 -> nx = (x + 1) % N;
                case 2 -> ny = (N + y - 1) % N;
                case 3 -> ny = (y + 1) % N;
            }
        }
        return new int[]{ny,nx};
    }
}
