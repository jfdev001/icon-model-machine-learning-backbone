# ICON Contribution Guidelines

## Introduction

ICON is simultaneously developed in several repositories. There is a _primary_ repository [icon](https://gitlab.dkrz.de/icon/icon) and several _secondary_ ones: [icon-nwp](https://gitlab.dkrz.de/icon/icon-nwp), [icon-mpim](https://gitlab.dkrz.de/icon/icon-mpim), etc. The default branch of the primary repository is also known as the release candidate (RC) of ICON.

## Communication

Use [issues](https://docs.gitlab.com/user/project/issues/) of the relevant ICON repository (fork) as a tool to communicate, track, and obtain information related to work items. Create an issue to report a bug you found, request a feature that needs to be implemented, or to discuss and coordinate with others on a particular topic. Issues can also be employed to address other work items, such as tracking the porting of a module to a different programming language or collecting information regarding an unexpected simulation result.

> **Note:** Prefer GitLab communication over private emails or messages to ensure information is searchable and accessible to a larger audience.

> **Note:** Avoid using [tasks](https://docs.gitlab.com/user/tasks/), as they add an extra level of hierarchy. To keep things simple, use issues only.

### Creating Issues

1. **Choose a template** (e.g. `bug-report`, `feature-request`) if available before writing an issue description.

2. Present issues clearly and include all relevant information so others can easily understand them.

3. Each issue should focus on **one actionable task** (e.g., one bug or feature request). If it becomes too lengthy or complex, break it into smaller issues and link them to the original one.

4. Always assign an issue to someone and communicate this clearly. Mention developers as needed, but keep it to a minimum.

5. Apply as many relevant [labels](https://docs.gitlab.com/user/project/labels) as possible. They help in classifying the issue (e.g. `Bug`, `Feature request`, `Discussion`) and provide additional context like machines used (e.g., `Levante`, `LUMI`), parts of ICON affected (e.g. `Ocean`, `CI`), or projects involved (e.g. `NextGEMS`, `WarmWorld`). Feel free to create new labels if the right ones don't exist.

> **Note:** Issues will be marked as `Stale` after two months of inactivity. Maintainers may close issues stalled for longer than two months. If the issue is still relevant, you're welcome to reopen it with a description following these guidelines.

## Code Contributions

### Coding style

We use [`pre-commit`](https://pre-commit.com) hooks to maintain a set of formatting and linting rules. Although there is a CI job that runs for each merge request and checks whether the contribution does not break the rules, we recommend registering the hooks in your local repository clone. This way, each commit undergoes the formatting and linking checks automatically.

We recommend installing `pre-commit` to a separate Python virtual environment using `pip`. For example, the following commands install the tool to the user's home directory:
```bash
python3 -m venv ~/pre-commit
~/pre-commit/bin/python3 -m pip install --upgrade pip
~/pre-commit/bin/python3 -m pip install pre-commit
```

You can now switch to the root of the repository and run the following command to register the hooks specified in [`.pre-commit-config.yaml`](/.pre-commit-config.yaml):
```bash
cd icon
~/pre-commit/bin/pre-commit install
```

From now on, each commit you make will be checked by a set of formatters and linters. Normally, the formatting tools are configured to modify the files in place. This means that if they fail, all you need to do is to accept the suggested changes and commit them:
```bash
git add .
git commit
```

Note that you will need to register the hooks for each fresh clone of the repository. Alternatively, you can follow [these instructions](https://pre-commit.com/#automatically-enabling-pre-commit-on-repositories) to configure `git` to register hooks automatically for each new clone of a repository that declares them.

### General Coding Rules

1. Avoid adding comments about future actions. For example,
    ```fortran
    ! Delete the subroutines below once the module is validated
    ```
    If necessary, include a reference to an issue that provides the progress status (e.g., an issue on the validation of the aforementioned module).
2. Do not add commented-out code to the codebase, as it produces maintenance and development overhead.

### Merge Requests

1. **Choose a template** (e.g. `nwp-feature`) if available before writing a merge request description.

2. Make the merge request title concise (titles become the first line of the commit message when the merge requests are accepted, together with a repository-specific prefix, e.g. `[mpim] `, `[nwp] `, etc.).

3. Please, adhere to the following recommendations for the merge requests **short** descriptions, which will become part of the commit message when the merge request is accepted:

    - use simple English in the active form (e.g. this implements A, updates B);

    - avoid special Markdown symbols and prefer plain ASCII, the message should read well in the terminal;

    - keep it short (excluding details, descriptions are appended to the merge request commit message);

    - do not reference issues and merge requests unless necessary (if referencing is necessary, make sure the reference contains the name and the namespace of the respective repository, e.g. `icon/icon#<issue-id>` and `icon/icon!<mr-id>`);

    - break the lines to make them no longer than 80 characters.

    > **Note:** The recommendations above apply to the **short** descriptions only. There are no restrictions for the **detailed** section.

4. The lists of co-authors in merge requests are generated automatically based on the authorship of the commits in the source branches. Please ensure that the commits in the source branch have the correct authorship with the correct email addresses (they can be [automatically-generated private commit emails](https://docs.gitlab.com/user/profile/#use-an-automatically-generated-private-commit-email)). If some commits have the wrong authorship, you can provide the list of co-authors using the following format:
    ```
    Co-authored-by: First-Name Second-Name <email.address@example.de>
    Co-authored-by: Another Name <another.address@example.com>
    ```

> **Note:** Tagging issues shows activity and prevents them from becoming `Stale`. Writing `Closes #<issue-ID>` in a merge request description automatically closes the issue when the request is merged.

### Documentation

The [ICON documentation](https://docs.icon-model.org) is automatically generated from the content of the subdirectory `doc/www/`. Please consider extending and updating the documentation with each merge request.
