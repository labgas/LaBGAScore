# conference_presentations — posters and talks, shared openly

LaBGAS posters and talks, one PDF each, shared with the community. Every file here is
also linked automatically from its first author's entry on the
[Team page](https://labgas.github.io/team/) of the lab website.

This is the one folder in LaBGAScore that holds no code.

## Naming convention

```
conference_presentations/<year>/<meeting>/<First_Author_Name>_<topic>_<meeting>_<year>.pdf
```

For example:

```
conference_presentations/2026/FNM/Maaike_Van_Den_Houte_TSPO_IBS_FNM_2026.pdf
```

- **`<year>`** — four digits, the year of the meeting.
- **`<meeting>`** — a short meeting name, such as `FNM`, `DDW` or `UEG_Week`. The folder name and
  the end of the file name must be the same. Use `_` instead of spaces.
- **`<First_Author_Name>`** — your name **exactly as it appears on the Team page**, with `_` for
  spaces. Accents and capital letters do not matter, so `Danielle_Jongen` and `Daniëlle_Jongen`
  both work. Hyphens may stay: `Anne-Sophie_Nobel`.
- **`<topic>`** — a few words, `_` between them. They become the link text, so keep it short:
  `TSPO_IBS` appears as **TSPO IBS · FNM 2026**.
- **PDF only.** Export slides to PDF before uploading.

## How it reaches the website

1. Commit and push the PDF to this repository. A file that only exists on your computer or the
   server is not picked up.
2. Once a day, the website looks at this folder and adds a **Conference presentations** link
   under the first author's name. Expect it within a day.

Only the first author is linked, because the file name names nobody else.

If your link does not appear after a day, the file name probably does not match. The
[website's daily run](https://github.com/labgas/labgas.github.io/actions) lists every file it could
not place, with the reason. The usual reasons are:
- the name differs from the one on the Team page;
- the meeting or year at the end of the file name differs from the folders;
- the first author is not on the Team page yet.

Fix the name and push again. To rename a file, use `git mv` so its history is kept.

## Licence

Unlike the code in this repository (GPL-3.0), the presentations in this folder are shared under
[Creative Commons Attribution 4.0 (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/):
anyone may reuse them with credit to the authors.

Before uploading, make sure every co-author agrees to sharing, and that the PDF contains nothing
you may not share. That includes unpublished data you want to keep back, third-party figures,
and anything that could identify a participant.
