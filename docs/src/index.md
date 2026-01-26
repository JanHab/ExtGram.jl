```@meta
CurrentModule = ExtGram
```

```@eval
using Markdown
Markdown.parse(read(joinpath(@__DIR__, "..", "..", "README.md"), String))
```
