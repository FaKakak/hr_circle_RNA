import java.util.Random;

public class Environment {
    public static final int SIDE = 30;
    public static final int TOTAL_MATERIAL = 80000;
    public static final int CELLNUM = SIDE*SIDE;
    int[][] raw_arr;
    RNA[][][] room_head;
    Random random;

    public Environment(Random random) {
        raw_arr = new int[SIDE][SIDE]; // every cell represents the quantity of raw material
        room_head = new RNA[2][SIDE][SIDE]; // every cell represents a RNA linkedlist
        this.random = random;
        // initialize the head of every RNA linkedlist
        for (int i = 0; i < 2; i++) {
            for(int j=0;j<SIDE;j++){
                for(int k=0;k<SIDE;k++){
                    room_head[i][j][k] = new RNA();
                }
            }
        }

        // initialize the raw_arr
        for(int i=0;i<TOTAL_MATERIAL;i++){
            int x = random.nextInt(SIDE);
            int y = random.nextInt(SIDE);
            raw_arr[y][x]++;
        }
    }
}
