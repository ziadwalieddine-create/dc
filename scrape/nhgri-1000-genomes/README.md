# NHGRI: 1000 Genomes Project (archived page scrape)

Source: https://www.genome.gov/27528684/1000-genomes-project
Page status: archived by NHGRI ("historical reference purposes only"); last updated 2017-10-04. Fetched 2026-10-25.

This is descriptive/programmatic material about the project. It contains no sequence or variant data.

## Files
- `page.html`: raw HTML of the project page
- `content.txt`: rendered main-content text (HTML comments removed)
- `1000genomes.json`: structured record: rationale, pilot projects, contacts, both tutorials (speaker, topic, affiliation, slide file, video URL where one exists), external links, and content hidden in HTML comments
- `linked_ashg-2010-...html/.txt`, `linked_ichg-2011-...html/.txt`: the two tutorial pages
- `linked_adam-felsenfeld-phd.html/.txt`: program director profile
- `assets/hh_1000gen.jpg`: hero image
- `assets/tutorials/ashg2010/` (6) and `assets/tutorials/ichg2011/` (5): slide PDFs, each with a `.txt` extracted via `pdftotext -layout` (slide text only; figures are not captured)

## Gaps and caveats
- Lisa Brooks's staff profile (`/staff/Lisa-Brooks-PhD`) returns HTTP 403, so it is not included. Her name, role and email are recorded from the project page.
- Videos are not downloaded. ASHG 2010 YouTube URLs are in the JSON. The ICHG 2011 page's "Watch" links are commented out in its source and reuse ASHG video IDs that do not match those speakers, so none are recorded for ICHG.
- The ICHG structural-variants slide link on the page uses the directory spelling `ICGH-...` (404). The same file at `ICHG-1000GenomesTutorial/Structural_Variation.pdf` works and is what was saved.
- An "Upcoming Meeting" block (July 2012 community analysis meeting) sits inside an HTML comment, so it is not shown on the live page; it is kept under `hidden_commented_out_content` in the JSON.
- External sites linked from the page (1000genomes.org, HapMap, YouTube) are listed but not scraped.
