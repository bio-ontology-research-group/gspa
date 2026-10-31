package gspa.predictor.genecalling

/** Prodigal prokaryotic gene caller, supporting single and metagenomic modes. */
class ProdigalCaller extends AbstractProdigalCaller {
    @Override
    String getName() { 'prodigal' }

    @Override
    protected String versionArgument() { '-v' }
}
