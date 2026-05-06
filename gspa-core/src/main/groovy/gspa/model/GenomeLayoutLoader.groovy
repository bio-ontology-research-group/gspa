package gspa.model

/**
 * Parse a layout TSV written by {@code benchmark/panel/make_layout.py}.
 *
 * <p>Columns (with header):
 * <pre>protein_id  contig  start  end  strand</pre></p>
 *
 * <p>Rows with missing / malformed columns are skipped with a warning.</p>
 */
class GenomeLayoutLoader {

    static GenomeLayout load(File tsv, String genomeId = null) {
        GenomeLayout layout = new GenomeLayout(genomeId)
        tsv.withReader { r ->
            String header = r.readLine()
            if (header == null) return
            // expected header; tolerate lowercase / extra whitespace
            int[] col = resolveColumns(header)
            String line
            while ((line = r.readLine()) != null) {
                if (!line || line.startsWith('#')) continue
                String[] parts = line.split('\t', -1)
                if (parts.length <= col[4]) continue
                try {
                    int start = Integer.parseInt(parts[col[2]].trim())
                    int end = Integer.parseInt(parts[col[3]].trim())
                    Strand strand = Strand.fromSymbol(parts[col[4]].trim())
                    layout.add(new ProteinLocus(
                        genomeId: genomeId,
                        proteinId: parts[col[0]].trim(),
                        contig: parts[col[1]].trim(),
                        start: start,
                        end: end,
                        strand: strand,
                    ))
                } catch (NumberFormatException ignored) {
                    // skip malformed row
                }
            }
        }
        layout.finishLoading()
        layout
    }

    private static int[] resolveColumns(String header) {
        List<String> names = header.split('\t').collect { it.trim().toLowerCase() }
        int[] out = new int[5]
        out[0] = names.indexOf('protein_id')
        out[1] = names.indexOf('contig')
        out[2] = names.indexOf('start')
        out[3] = names.indexOf('end')
        out[4] = names.indexOf('strand')
        for (int i = 0; i < 5; i++) if (out[i] < 0) out[i] = i
        out
    }
}
