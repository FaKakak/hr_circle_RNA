public class c2_frag {
    int start;
    int length;
    c2_frag next;
    c2_frag prior;

    public c2_frag() {
        this.start = 0;
        this.length = 0;
        this.next = this;
        this.prior = this;
    }
    public c2_frag(int start,int length) {
        this.start = start;
        this.length = length;
        this.next = this;
        this.prior = this;
    }
    public void addAfter(c2_frag newFrag){
        newFrag.next = next;
        newFrag.prior = this;
        next.prior = newFrag;
        next = newFrag;
    }
}
