<!--
# YAVSAP pull request

Many thanks for contributing to YAVSAP!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the develop branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/ksumngs/yavsap/tree/master/.github/CONTRIBUTING.md)
-->
<!-- markdownlint-disable ul-indent -->
<!-- markdownlint-disable line-length -->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
  - [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/ksumngs/yavsap/tree/master/.github/CONTRIBUTING.md)
  - [ ] If necessary, also make a PR on the yavsap _branch_ on the [ksumngs/nf-test-datasets](https://github.com/ksumngs/nf-test-datasets) repository.
- [ ] Make sure your code lints (`nf-core lint`).
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
