package gspa.predictor.genecalling

/** Pyrodigal prokaryotic gene caller. Configuration may select Prodigal if unavailable. */
class PyrodigalCaller extends AbstractProdigalCaller {
    @Override
    String getName() { 'pyrodigal' }

    @Override
    protected String versionArgument() { '-V' }
}
