---
description: "BEACH の GitHub Issue 草案を作る"
argument-hint: "<problem description/log/config>"
---

`beach-issue-report` agent を使い、問題内容を GitHub Issue draft に整理してください。

出力には次を含めてください:
- title
- type: Bug report / Improvement request / Feature request / Documentation request
- body
- suggested labels
- missing information

token、private path、巨大 CSV 全文は貼らず、短い excerpt と file listing を優先してください。
