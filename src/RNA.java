import java.lang.System;
import java.util.Arrays;

public class RNA {
    // 跟RNA有关的宏定义都写在这儿了
    public static final int MAX_RNA_LENGTH = 100;
    public static final int MAX_CHAR_LENGTH=20;
    public static final char[] NRSEQ = {'U','G','A','U','G','C','A','G'};
    public static final char[] CTR1SEQ = {'A','C','U','G','A','C','G','U'};
    public static final char[] CTR2SEQ = {'A','A','C','G','C','U','C','G'};
    public static final char[] CTR3SEQ = {'C','G','A','U','C','A','A','U'};
    public static final char[] HRSEQ = {'A','G','U','C'};
    public static final char[] INOCUSEQ = NRSEQ;
    public static final char[] INOCUSEQ1 = CTR1SEQ;
    public static final char[] INOCUSEQ2 = CTR2SEQ;
    public static final char[] INOCUSEQ3 = CTR3SEQ;


    char[][] information;
    int length1;
    int length2;
    int type1;
    int type2;

    RNA next;
    RNA prior;
    c2_frag chain2;

    public RNA() {
        information = new char[2][MAX_RNA_LENGTH];
        length1=0;
        length2=0;
        type1=0;
        type2=0;
        next = this;
        prior = null;
        chain2= new c2_frag();
        for (int i = 0; i < MAX_RNA_LENGTH; i++) {
            information[0][i] = '0';
            information[1][i]='0';
        }

    }
    public RNA(char[] bond1,int type) {
        information = new char[2][MAX_RNA_LENGTH];
        length1= bond1.length;
        length2=0;
        type1=type;
        type2=0;
        next = this;
        prior = null;
        chain2= new c2_frag();
        for (int i = 0; i < MAX_RNA_LENGTH; i++) {
            information[0][i] = '0';
            information[1][i]='0';
        }
        System.arraycopy(bond1, 0, information[0], 0, bond1.length);
    }

    public void addAfter(RNA newRNA){// 把newRNA添加到该RNA后边
        RNA oldNext = next;
        next = newRNA;
        newRNA.prior = this;
        newRNA.next = oldNext;
        if(oldNext.prior!=null) oldNext.prior = newRNA;
    }
    public void removeThis(){
        // 游离的RNA不用动
        if(prior==null) {
            if(next == this) return ;
        }
        if(next.prior!=null) {
            // 链表头下边不止一个RNA
            next.prior = prior;
        }
        prior.next = next;
        next = this;
        prior = null;
    }
}
