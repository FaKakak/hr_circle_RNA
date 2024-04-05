import java.util.Random;

public class Inoculate_generator {
    public static final int INOCUNUM = 20;
    Environment environment;
    int[][] raw_arr;
    RNA[][][] cell_head;
    Random random;

    public Inoculate_generator(Environment environment,Random random) {
        this.random = random;
        this.environment = environment;
        this.raw_arr = environment.raw_arr;
        this.cell_head = environment.room_head;
    }

    public void inoculate(int h,int type){
        for (int i = 0; i < INOCUNUM; i++) {
            inoculate_helper(h,RNA.INOCUSEQ,type);
            inoculate_helper(h,RNA.INOCUSEQ1,type);
            inoculate_helper(h,RNA.INOCUSEQ2,type);
            inoculate_helper(h,RNA.INOCUSEQ3,type);
        }
    }

    private void inoculate_helper(int h,char[] seq,int type){
        int x = random.nextInt(Environment.SIDE);
        int y = random.nextInt(Environment.SIDE);
        RNA newRNA = new RNA(seq,type);
        cell_head[h][y][x].addAfter(newRNA);
    }
}
