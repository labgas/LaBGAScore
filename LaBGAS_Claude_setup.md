# LaBGAS Claude setup

How to go from receiving Lukas's e-mail invitation to our Claude **Team** account to being
fully operational, on both your own Windows machine (Claude desktop app) and the shared
LaBGAS Linux server (Claude Code CLI), where everyone logs in with their own KU Leuven
account and home directory.

## Contents

- [1. Accept the Team invitation](#1-accept-the-team-invitation)
- [2. Windows desktop app](#2-windows-desktop-app)
- [3. Claude Code CLI on the shared Linux server](#3-claude-code-cli-on-the-shared-linux-server)
- [Getting help](#getting-help)

## 1. Accept the Team invitation

1. Lukas sends a Team invite to your KU Leuven e-mail address from the Claude admin console.
2. Open the invite e-mail and follow the link. If you don't already have a personal Claude
   account under that address, create one; otherwise sign in with your existing one.
3. Accepting the invite joins you to the LaBGAS organization/Team. This is the **same
   account** you use everywhere below — desktop app and CLI both authenticate against it,
   there is nothing separate to set up per tool.

## 2. Windows desktop app

1. **Install.** Download the Claude desktop app from [claude.ai/download](https://claude.ai/download)
   and install it as usual.
2. **Sign in.** Use the KU Leuven e-mail address you accepted the Team invite with. You
   should land in the LaBGAS organization automatically.
3. **Connect Slack.** In the app's connector/integration settings, add the Slack connector
   and authorize it against the LaBGAS Slack workspace. This is what lets Claude read/search
   our channels and draft messages when you ask it to.
4. **Connect GitHub.** Add the GitHub connector and authorize it for the `labgas` GitHub
   organization (at minimum, `LaBGAScore` and `CANlab_help_examples`, plus any `proj_xxx`
   repos you work on). This is what lets Claude open PRs, read issues, and check CI for you.
5. **MATLAB Agentic Toolkit (optional, recommended if you run MATLAB locally on Windows).**
   Already fully documented — don't redo it here. See step 7 of
   ["Before you start"](LaBGAS_fMRI_analysis_workflow.md#before-you-start) in the workflow
   document for the install commands and the LaBGAS-recommended skill-group subset.

You're operational once you can open a chat, see the LaBGAS organization in the app, and
(if you set them up) see Slack/GitHub listed as connected.

## 3. Claude Code CLI on the shared Linux server

This is a **separate install per person**, even though you're all on the same machine.
Nothing here needs root or a linuxteam ticket — each of you installs into your own home
directory, and your login, settings, and conversation memory all live there too, isolated
from everyone else logged into the same server.

> **Why per-user, not system-wide:** a shared, root-owned install would need linuxteam to
> set up and to keep updated (Claude Code updates itself in place, which a non-root user
> can't do to a root-owned binary), for a saving of one copy-paste command per person.
> Per-user costs nothing but that one command, with no admin dependency.

1. **Connect to the server** via X2go (or SSH) with your own KU Leuven institutional login,
   as usual.

2. **Install the CLI** into your own home directory:

   ```bash
   curl -fsSL https://claude.ai/install.sh | bash
   ```

   This installs a native binary under `~/.local/share/claude/versions/` and symlinks it
   from `~/.local/bin/claude` — no Node.js/npm required.

3. **Make sure `~/.local/bin` is on your `PATH`.** Most shells already include it; check
   with:

   ```bash
   which claude
   ```

   If that prints nothing, add `export PATH="$HOME/.local/bin:$PATH"` to your `~/.bashrc`
   (or `~/.zshrc`), then open a new shell.

4. **Log in.** Run:

   ```bash
   claude
   ```

   On first run it walks you through a one-time login: it prints a URL, which you open in
   a browser **on your own computer** (the server itself has no GUI browser in a plain
   SSH/headless session), sign in with the same KU Leuven Team account from step 1, and
   authorize. Return to the terminal and you're logged in. Your credentials are stored
   under your own `~/.claude*`, so this is entirely independent of whoever else is logged
   into the server at the same time.

5. **Trust the repo the first time you open one.** The first time you run `claude` inside
   a new project directory (e.g. `/data/master_github_repos/LaBGAScore` or a `proj_xxx`
   dataset's `code` subdataset), it asks whether you trust the files there. Say yes for
   LaBGAS repos — this is also when Claude picks up that repo's `CLAUDE.md`/`README.md`
   conventions automatically.

6. **Optional: MATLAB Agentic Toolkit here too.** If you run MATLAB from this server (see
   the headless-run instructions in
   ["Before you start"](LaBGAS_fMRI_analysis_workflow.md#before-you-start) of the workflow
   document), the same MCP server and skill groups documented there work on Linux as well
   — the install command and recommended skill-group subset are identical; only the
   `setupAgenticToolkit("install")` step needs to be run from MATLAB on the server rather
   than on Windows.

You're operational once `which claude` resolves, `claude` starts without prompting you to
log in again, and you can open a chat from inside a LaBGAS repo directory.

## Getting help

- CLI/desktop app usage questions: ask Claude itself (`/help` in the CLI), or see the
  [Claude Code documentation](https://docs.claude.com/en/docs/claude-code).
- Team account / seat / billing issues: contact Lukas.
- Server access issues (X2go, your KU Leuven login not working, needing root for something
  unrelated to Claude): [linuxteam.gbiomed@kuleuven.be](mailto:linuxteam.gbiomed@kuleuven.be)
  (or, better, an ICTS helpdesk ticket with Lukas in cc) — see step 1 of
  ["Before you start"](LaBGAS_fMRI_analysis_workflow.md#before-you-start) in the workflow
  document.
