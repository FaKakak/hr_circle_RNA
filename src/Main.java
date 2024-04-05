import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import static javax.print.attribute.standard.MediaTray.SIDE;

public class Main {
    public static final int STEPNUM = 1500000;
    public static final int INOCUSTEP = 10000;
    public static final int STAREC = 0;
    public static final int RECINT = 10000;
    public static final int SEED = 825;

    public static void main(String[] args) {
        Random random = new Random(SEED);
        long startTime = System.currentTimeMillis();
        int h=0;
        Environment e = new Environment(random);
        Unit_case_generator unitCaseGenerator = new Unit_case_generator(e,random);
        Record_generator recordGenerator = new Record_generator(e);
        Inoculate_generator inoculateGenerator = new Inoculate_generator(e,random);
        for (int i = 0; i <= STEPNUM; i++) {
            if(i==INOCUSTEP) {
                inoculateGenerator.inoculate(h,0);
            }
            if(i>=STAREC&&i%RECINT==0){
                recordGenerator.record(h,i);
            }
            unitCaseGenerator.unit_case(h);
            h=h^1;
        }
        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
        try (PrintWriter writer = new PrintWriter(new FileWriter("Java.txt", true))) {
            writer.printf("程序运行时间："+duration*0.001+" s, SD:"+SEED);
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        }
    }
}