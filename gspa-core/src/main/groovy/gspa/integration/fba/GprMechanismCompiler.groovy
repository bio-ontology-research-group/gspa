package gspa.integration.fba

import groovy.transform.CompileStatic

/**
 * Stage A1: Compiles reaction Gene-Protein-Reaction (GPR) boolean expressions
 * and multi-subunit EC mappings into Disjunctive Normal Form (DNF) mechanism clauses
 * M(j) = { k_1, ..., k_m }, where each mechanism k = { r_1, ..., r_L } is a conjunction
 * (AND) of required catalytic/subunit roles.
 */
@CompileStatic
class GprMechanismCompiler {

    static class MechanismClause {
        final List<String> requiredRoles
        final boolean isHeteromericComplex

        MechanismClause(List<String> requiredRoles, boolean isHeteromericComplex) {
            this.requiredRoles = requiredRoles
            this.isHeteromericComplex = isHeteromericComplex
        }
    }

    /**
     * Compile a reaction's mapped EC numbers and optional raw GPR string into DNF clauses.
     * Supports nested boolean GPR expressions with parentheses, 'and', and 'or'
     * (e.g. "(b0001 and b0002) or b0003"), expanding them into Disjunctive Normal Form (DNF).
     *
     * @param mappedEcs List of complete 4-digit EC numbers mapped to reaction j
     * @param rawGpr    Raw boolean GPR string (e.g. "(b0001 and b0002) or b0003")
     * @return List of MechanismClause representing alternative isozyme/complex mechanisms
     */
    static List<MechanismClause> compileMechanisms(List<String> mappedEcs, String rawGpr = '') {
        List<String> cleanEcs = new ArrayList<String>()
        if (mappedEcs != null) {
            for (String ec : mappedEcs) {
                if (ec == null) continue
                String c = ec.replace('EC:', '').trim()
                if (!c.isEmpty() && !c.contains('-')) {
                    cleanEcs.add(c)
                }
            }
        }
        String gprStr = (rawGpr ?: '').trim()
        if (!gprStr.isEmpty()) {
            List<List<String>> dnf = parseGprToDnf(gprStr)
            if (!dnf.isEmpty()) {
                Set<String> cleanEcSet = new LinkedHashSet<String>(cleanEcs)
                boolean tokensMatchRoles = cleanEcSet.isEmpty()
                if (!tokensMatchRoles) {
                    for (List<String> c : dnf) {
                        for (String tok : c) {
                            String normTok = tok.replace('EC:', '').trim()
                            if (cleanEcSet.contains(normTok)) {
                                tokensMatchRoles = true
                                break
                            }
                        }
                        if (tokensMatchRoles) break
                    }
                }
                if (tokensMatchRoles || dnf.size() > 1) {
                    List<MechanismClause> gprClauses = new ArrayList<MechanismClause>(dnf.size())
                    for (List<String> c : dnf) {
                        List<String> normClause = new ArrayList<String>(c.size())
                        for (String tok : c) {
                            normClause.add(tok.replace('EC:', '').trim())
                        }
                        gprClauses.add(new MechanismClause(normClause, normClause.size() > 1))
                    }
                    return gprClauses
                }
                if (dnf.get(0).size() > 1 && cleanEcs.size() > 1) {
                    return Collections.singletonList(new MechanismClause(cleanEcs, true))
                }
            }
        }
        if (cleanEcs.isEmpty()) {
            return Collections.emptyList()
        }

        List<MechanismClause> clauses = new ArrayList<MechanismClause>(cleanEcs.size())
        for (String ec : cleanEcs) {
            clauses.add(new MechanismClause(Collections.singletonList(ec), false))
        }
        return clauses
    }

    /**
     * Parse a boolean GPR expression with parentheses, AND, and OR into Disjunctive Normal Form
     * (list of conjunction clauses, where each clause is a list of required gene/role tokens).
     */
    static List<List<String>> parseGprToDnf(String expr) {
        if (expr == null || expr.trim().isEmpty()) return Collections.emptyList()
        String normalized = expr.replace('(', ' ( ').replace(')', ' ) ')
        String[] rawTokens = normalized.trim().split('\\s+')
        List<String> tokens = new ArrayList<String>(rawTokens.length)
        for (String t : rawTokens) {
            if (!t.isEmpty()) tokens.add(t)
        }
        if (tokens.isEmpty()) return Collections.emptyList()
        int[] pos = [0] as int[]
        return parseOrExpr(tokens, pos)
    }

    private static List<List<String>> parseOrExpr(List<String> tokens, int[] pos) {
        List<List<String>> left = parseAndExpr(tokens, pos)
        while (pos[0] < tokens.size() && tokens.get(pos[0]).equalsIgnoreCase('or')) {
            pos[0]++
            List<List<String>> right = parseAndExpr(tokens, pos)
            List<List<String>> combined = new ArrayList<List<String>>(left.size() + right.size())
            combined.addAll(left)
            combined.addAll(right)
            left = combined
        }
        return left
    }

    private static List<List<String>> parseAndExpr(List<String> tokens, int[] pos) {
        List<List<String>> left = parsePrimary(tokens, pos)
        while (pos[0] < tokens.size() && tokens.get(pos[0]).equalsIgnoreCase('and')) {
            pos[0]++
            List<List<String>> right = parsePrimary(tokens, pos)
            left = crossProductAnd(left, right)
        }
        return left
    }

    private static List<List<String>> parsePrimary(List<String> tokens, int[] pos) {
        if (pos[0] >= tokens.size()) return Collections.emptyList()
        String tok = tokens.get(pos[0])
        if (tok == '(') {
            pos[0]++
            List<List<String>> inner = parseOrExpr(tokens, pos)
            if (pos[0] < tokens.size() && tokens.get(pos[0]) == ')') {
                pos[0]++
            }
            return inner
        }
        pos[0]++
        return Collections.singletonList(Collections.singletonList(tok))
    }

    private static List<List<String>> crossProductAnd(List<List<String>> a, List<List<String>> b) {
        if (a.isEmpty()) return b
        if (b.isEmpty()) return a
        List<List<String>> out = new ArrayList<List<String>>(a.size() * b.size())
        for (List<String> ca : a) {
            for (List<String> cb : b) {
                Set<String> merged = new LinkedHashSet<String>(ca)
                merged.addAll(cb)
                out.add(new ArrayList<String>(merged))
            }
        }
        return out
    }
}
