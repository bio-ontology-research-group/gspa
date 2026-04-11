package gspa.model

enum Strand {
    PLUS('+'),
    MINUS('-'),
    NONE('.')

    final String symbol

    Strand(String symbol) {
        this.symbol = symbol
    }

    static Strand fromSymbol(String s) {
        switch (s) {
            case '+': return PLUS
            case '-': return MINUS
            default: return NONE
        }
    }

    boolean isSameStrand(Strand other) {
        this == other && this != NONE
    }
}
