# ROTIFER — dossiê de código legado candidato a remoção

**Repositório:** `leepusp/rotifer` (clone completo de https://github.com/leepusp/rotifer.git)  
**HEAD analisado:** `522f2871725f` — 2026-09-09 — "hot fix"  
**Commits no histórico:** 1331  
**Data da auditoria:** 2026-09-11  
**Natureza:** READ-ONLY. Nenhum arquivo foi editado, movido ou removido; nenhum branch, commit, issue ou PR foi criado.

**Total de candidatos: 230** (Alta: 97 · Média: 72 · Baixa/investigar: 61)

---

## 0. Escopo, caminhos e método

### Caminho de `devel/alpha`

Resolvido a partir da raiz do clone **antes** da varredura:

```
lib/rotifer/devel/alpha/            <- sandbox pessoal dos desenvolvedores
lib/rotifer/devel/alpha/genome/
lib/rotifer/devel/alpha/sequence/
```

Ou seja, `devel/alpha` vive **dentro** da árvore do núcleo (`lib/rotifer/`). Toda contagem de call site em "núcleo" exclui explicitamente o prefixo `lib/rotifer/devel/alpha/`. `lib/rotifer/devel/beta/` **não** é sandbox e foi contado como núcleo.

### Arquivos enumerados

| Categoria                                                              | Arquivos |
| ---------------------------------------------------------------------- | -------- |
| Python (`.py`, `.ipynb`, shebang `python*`, cópias `.bkp`/`.quebrada`) | 199      |
| Perl (`.pl`, `.pm`, shebang `perl`)                                    | 164      |
| Shell (`.sh`, `.zsh`, `.lib`, shebang `bash`/`sh`/`zsh`)               | 58       |
| **Subtotal código**                                                    | **421**  |
| Outros arquivos rastreados (dados, docs, configs, binários)            | 171      |
| **Total rastreado por `git ls-files`**                                 | **592**  |

| Definições extraídas               | Quantidade |
| ---------------------------------- | ---------- |
| Python — `def` (funções)           | 1115       |
| Python — `def` (métodos de classe) | 753        |
| Python — `class`                   | 148        |
| Perl — `sub`                       | 547        |
| Shell — funções                    | 59         |
| **Total**                          | **2622**   |

Distribuídas em **323 arquivos**. Extração Python por `ast` (`ast.parse` + visitor, com `lineno`/`end_lineno` reais); Perl por varredura de `^\s*sub\s+NOME` com resolução de `package` e casamento de chaves para o fim do corpo; Shell por `nome() {` / `function nome {` com o mesmo casamento de chaves.

### Arquivos que não compilam (evidência dura de abandono)

Cinco arquivos Python **não passam por `ast.parse`** — isto é, não podem ser importados nem executados em nenhuma versão de Python 3:

| Arquivo                                                      | Erro                                |
| ------------------------------------------------------------ | ----------------------------------- |
| `lib/rotifer/db/neighbors.py`                                | f-string: expecting '}' (linha 131) |
| `lib/rotifer/io/base.py`                                     | '(' was never closed (linha 6)      |
| `share/rotifer/snakemake/rps/Snakefile.bkp`                  | invalid syntax (linha 6)            |
| `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp2`        | invalid syntax (linha 5)            |
| `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp_working` | invalid syntax (linha 5)            |

Dois deles estão no **núcleo** (`lib/rotifer/io/base.py`, `lib/rotifer/db/neighbors.py`). Verificação manual em `lib/rotifer/io/base.py:6`: a chamada `super().__init__(` é aberta e nunca fechada. Em `lib/rotifer/db/neighbors.py:131`: f-string `{cursor.fetchone()[0]` sem `}` de fechamento. Isto significa que **todo símbolo definido nesses dois arquivos é inalcançável em tempo de execução** — nenhum `import rotifer.io.base` pode ter sucesso.

### Como cada eixo foi medido

**Eixo 1 — Idade.** `git blame --line-porcelain -w` em cada arquivo com definições (o `-w` já descarta reindentação/espaço em branco). Para o intervalo de linhas de cada função coletou-se: commit mais recente, commit mais antigo sobrevivente, nº de commits distintos e nº de autores distintos. A data de **introdução** foi confirmada independentemente por pickaxe (`git log -S "def nome(" --reverse`). Quando os dois métodos divergem, ambos aparecem na ficha.

**Commits de ruído identificados e descartados do cálculo de "última alteração significativa":**

| Commit    | Data       | Natureza                                                                    | Por que é ruído                                                                           |
| --------- | ---------- | --------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `b0e3490` | 2019-05-02 | Modified rotifer programs name                                              | 26 renames puros em `bin/` (`acc2operon.py` → `acc2operon`), 0 linhas de lógica alteradas |
| `462710b` | 2019-05-02 | Removed pycache                                                             | 80 arquivos, 0 inserções / 0 deleções — só `.pyc`                                         |
| `8303356` | 2019-05-02 | Revert "Removed pycache"                                                    | 80 arquivos, 0 inserções / 0 deleções — reversão do anterior                              |
| `29bd8c1` | 2019-03-15 | Modified directory structure                                                | 168 arquivos, +14510/−21 — movimentação de diretórios registrada como adição              |
| `3364c1a` | 2019-03-15 | Fixing sys.path in scripts. Removed duplicated files.                       | 156 arquivos, +28/−13789 — limpeza pós-reestruturação                                     |
| `81f751c` | 2019-05-02 | Modified structure of rotifer module                                        | 50 arquivos — reestruturação do pacote                                                    |
| `8214b5d` | 2019-09-18 | First attempt at removing hard-coded references to Kaihami's home directory | 21 arquivos, +14/−12 — substituição mecânica de caminhos                                  |

**Commits de importação em massa** (marcam a _introdução_ do código, não uma alteração de conteúdo) — também descartados do cálculo de "última alteração significativa":

| Commit    | Data       | Natureza                                                                                                 |
| --------- | ---------- | -------------------------------------------------------------------------------------------------------- |
| `3adf523` | 2019-03-08 | "First commit adding all files, including the ones that I should remove" — 397 arquivos, +632.647 linhas |
| `0602ce1` | 2019-05-02 | "Merging old rotifer Perl code and scripts." — 249 arquivos, +72.161 linhas (toda a árvore `perl/`)      |

Os commits de licenciamento (`3863a74`, `834b292`, `765a61b`, 2026) foram inspecionados e **não tocam nenhum arquivo de código** — só `LICENSE.md`, `CODE_OF_CONDUCT.md`, `.github/` e `doc/licenses/`. Não houve, portanto, ruído de licença a descartar.

**Eixo 2 — Uso.** Índice de identificadores construído lendo **os 586 arquivos de texto rastreados**, linha a linha, e registrando toda ocorrência de cada identificador com `arquivo:linha`. O arquivo `tags` (índice ctags gerado, 341 KB, versionado por engano) foi **excluído do índice** para não inflar contagens. As ocorrências são particionadas em quatro baldes, **nunca somados**:

- **(a) núcleo** — `lib/` exceto `lib/rotifer/devel/alpha/`
- **(b) CLI** — `bin/` e `perl/`
- **(c) código pessoal** — `lib/rotifer/devel/alpha/`
- **(d) testes e documentação** — `test/` e `doc/`

(Um quinto balde, `outros` — `etc/`, `share/`, `sbin/`, raiz —, é registrado nas fichas individuais mas não entra na coluna de quatro casas da tabela.)

A linha de definição e o corpo da própria função são descontados. **Qualquer função com contagem > 0 em (a) ou (b) foi removida da lista de candidatos, independentemente da idade** — foram 1.014 definições eliminadas por esse critério.

**Limite conhecido deste método:** a correspondência é por nome simples, não por resolução de escopo. Isso torna a medida _conservadora na direção segura_ — homônimos inflam a contagem de uso e tiram funções da lista, nunca o contrário. Onde um nome tem mais de uma definição no repositório, a ficha traz `NOME AMBÍGUO` e a confiança é rebaixada.

**Eixo 3 — Obsolescência.** Varredura do corpo de cada função (mais as 3 linhas anteriores, para pegar comentários de cabeçalho) por: marcadores textuais (`deprecated`, `obsolete`, `unused`, `legacy`, `TODO`, `FIXME`, `XXX`, `HACK`, `BROKEN`, `quebrad*`); idiomas Python 2 (`print` statement, `.iteritems/.itervalues/.iterkeys`, `.has_key`, `urllib2`, `optparse`, `import commands`, `string.join`, `except X, e`); caminhos absolutos hardcoded; blocos grandes de código comentado; e pertencimento a arquivo de backup versionado ou a arquivo que não compila.

### Exclusões obrigatórias aplicadas (marcadas NÃO REMOVER)

| Regra                                                                                  | Definições excluídas |
| -------------------------------------------------------------------------------------- | -------------------- |
| Símbolo exportado, importado ou definido em algum `__init__.py`                        | 192                  |
| Identificador citado em `doc/`, `README.md`, `CONTRIBUTING.md` ou `CODE_OF_CONDUCT.md` | 302                  |
| Entry point de ferramenta CLI (`main` em `bin/`/`sbin/`)                               | 0                    |
| **Total**                                                                              | **494**              |

A exclusão por documentação é deliberadamente generosa: basta o identificador aparecer em qualquer arquivo sob `doc/` (inclusive em prosa) para a função sair da lista. Isso custa alguns candidatos legítimos, mas respeita a restrição de compatibilidade retroativa do projeto.

### Funil completo — de 2.622 definições a 230 candidatos

| Etapa                                                                                            | Definições | Restam  |
| ------------------------------------------------------------------------------------------------ | ---------- | ------- |
| Total de definições extraídas                                                                    | —          | 2622    |
| −171 — Métodos `__dunder__` (invocados implicitamente pelo interpretador)                        | 171        | 2451    |
| −413 — Definidas em `test/` (funções `pytest`, coletadas por convenção, nunca chamadas por nome) | 413        | 2038    |
| −1 — Definidas em `doc/` (exemplos e templates)                                                  | 1          | 2037    |
| −494 — Exclusões obrigatórias (`__init__.py` / `doc/` / entry point CLI)                         | 494        | 1543    |
| −1014 — Com call site em núcleo (a) ou CLI (b) — **não candidatas por regra**                    | 1014       | 529     |
| −262 — Sem uso em (a)/(b) mas alteradas nos últimos 2 anos e sem nenhum sinal de obsolescência   | 262        | 267     |
| −37 — Falsos positivos descartados por análise manual (ver §5)                                   | 37         | 230     |
| **Candidatos finais**                                                                            | —          | **230** |

---

## 1. Alta confiança — remoção segura — 97 funções

## obs: TUDO NESSA LISTA (Alta confiança) FOI DELETADO.

Evidência convergente nos três eixos: sem call site em núcleo, CLI, testes ou documentação; histórico parado ou arquivo notoriamente descartado; e pelo menos um sinal objetivo de obsolescência. Nome não ambíguo, ou ambíguo apenas entre um arquivo e suas próprias cópias de backup.

#### 1. `TextIO` — Alta

- **Assinatura:** `TextIO` — python, class
- **Local:** `lib/rotifer/io/base.py:4`–4 (0 linhas)
- **Última alteração (bruta):** 2020-05-29 · **excluindo ruído e importação em massa:** 2020-05-29
- **Commits que tocaram o corpo (1):** `e25d017`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "class TextIO("`):** `e25d017` — 2020-05-29 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - ARQUIVO NÃO PARSEIA: '(' was never closed (linha 6)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo: o arquivo não compila, logo nenhum import dele jamais teve sucesso; nenhum workflow externo pode depender deste símbolo.

#### 2. `apply_selected_scheme(selected, color_data)` — Alta

- **Assinatura:** `apply_selected_scheme(selected, color_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:835`–839 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def apply_selected_scheme("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 3. `update_communities(n_clicks, serialized_G, color_data, selected_groups, min_connections, function_data)` — Alta

- **Assinatura:** `update_communities(n_clicks, serialized_G, color_data, selected_groups, min_connections, function_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:853`–883 (31 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_communities("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 4. `update_graph_store(min_connections, selected_groups, original_data, function_data)` — Alta

- **Assinatura:** `update_graph_store(min_connections, selected_groups, original_data, function_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:892`–909 (18 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_graph_store("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 5. `autoopen_seqio` — Alta

- **Assinatura:** `autoopen_seqio` (em `action`) — python, class
- **Local:** `lib/rotifer/seq/cli.py:33`–53 (21 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "class autoopen_seqio("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: todo
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 6. `_multiple_find_acc(sub_fi, gbk, db, send_end)` — Alta

- **Assinatura:** `_multiple_find_acc(sub_fi, gbk, db, send_end)` — python, function
- **Local:** `bin/acc2operon:624`–630 (7 linhas)
- **Última alteração (bruta):** 2019-03-08 · **excluindo ruído e importação em massa:** 2019-03-08
- **Commits que tocaram o corpo (1):** `3adf523`
- **Autores distintos (1):** kaihami
- **Commits de importação em massa no intervalo (descartados):** `3adf523`
- **Introdução (pickaxe `git log -S "def _multiple_find_acc("`):** `3adf523` — 2019-03-08 — kaihami · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 7. `exec_function(df, s, verbose)` — Alta

- **Decision:** remove file
- **Assinatura:** `exec_function(df, s, verbose)` — python, function
- **Local:** `bin/rnexplorer:201`–212 (12 linhas)
- **Última alteração (bruta):** 2019-03-08 · **excluindo ruído e importação em massa:** 2019-03-08
- **Commits que tocaram o corpo (1):** `3adf523`
- **Autores distintos (1):** kaihami
- **Commits de importação em massa no intervalo (descartados):** `3adf523`
- **Introdução (pickaxe `git log -S "def exec_function("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 8. `blockPrint()` — Alta

- **Assinatura:** `blockPrint()` — python, function
- **Local:** `bin/rotifer:31`–33 (3 linhas)
- **Última alteração (bruta):** 2019-05-06 · **excluindo ruído e importação em massa:** 2019-05-06
- **Commits que tocaram o corpo (1):** `5ff88be`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def blockPrint("`):** `81f751c` — 2019-05-02 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 9. `enablePrint()` — Alta

- **Assinatura:** `enablePrint()` — python, function
- **Local:** `bin/rotifer:35`–36 (2 linhas)
- **Última alteração (bruta):** 2019-05-06 · **excluindo ruído e importação em massa:** 2019-05-06
- **Commits que tocaram o corpo (1):** `5ff88be`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def enablePrint("`):** `81f751c` — 2019-05-02 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 10. `sub old_ncbi_gi_parser` — Alta

- **Assinatura:** `sub old_ncbi_gi_parser` — perl, sub
- **Local:** `bin/seqstats:69`–72 (4 linhas)
- **Última alteração (bruta):** 2020-01-13 · **excluindo ruído e importação em massa:** 2020-01-13
- **Commits que tocaram o corpo (2):** `0602ce1`, `f05b9ba`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "sub old_ncbi_gi_parser"`):** `f05b9ba` — 2020-01-13 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 11. `sub map_labels` — Alta

- **Decision:** reverse
- **Assinatura:** `sub map_labels` — perl, sub
- **Local:** `bin/treeutil:159`–173 (15 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "sub map_labels"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 12. `sub load_and_process` — Alta

- **Assinatura:** `sub load_and_process` — perl, sub
- **Local:** `bin/tupdate:21`–91 (71 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "sub load_and_process"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 13. `cluster_neighbor` — Alta

- **Assinatura:** `cluster_neighbor` — python, class
- **Local:** `lib/rotifer/cluster/cluster.py:14`–93 (80 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe):** [sem evidência]
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 14. `build_matrix(self)` — Alta

- **Assinatura:** `build_matrix(self)` (em `cluster_neighbor`) — python, method
- **Local:** `lib/rotifer/cluster/cluster.py:74`–93 (20 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def build_matrix("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 15. `add_arguments` — Alta

- **Assinatura:** `add_arguments` — python, class
- **Local:** `lib/rotifer/core/dev/cli2.py:11`–34 (24 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe):** [sem evidência]
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 16. `fetch_seq(seqs)` — Alta

- **Assinatura:** `fetch_seq(seqs)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:252`–274 (23 linhas)
- **Última alteração (bruta):** 2022-10-18 · **excluindo ruído e importação em massa:** 2022-10-18
- **Commits que tocaram o corpo (5):** `70879c7`, `db472c1`, `d24dddc`, `97bea48`, `658aff6`
- **Autores distintos (2):** Gianlucca Gonçalves Nicastro, Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def fetch_seq("`):** `70879c7` — 2021-12-10 — Gianlucca Gonçalves Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 17. `trim_unk_neigh(df, ann)` — Alta

- **Assinatura:** `trim_unk_neigh(df, ann)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:946`–959 (14 linhas)
- **Última alteração (bruta):** 2023-07-17 · **excluindo ruído e importação em massa:** 2023-07-17
- **Commits que tocaram o corpo (1):** `8da21f0`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def trim_unk_neigh("`):** `8da21f0` — 2023-07-17 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: xxx
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 18. `from_yaml()` — Alta

- **Assinatura:** `from_yaml()` (em `configuration`) — python, method
- **Local:** `lib/rotifer/devel/alpha/mapper.py:29`–30 (2 linhas)
- **Última alteração (bruta):** 2022-06-22 · **excluindo ruído e importação em massa:** 2022-06-22
- **Commits que tocaram o corpo (1):** `8c6a156`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def from_yaml("`):** `8c6a156` — 2022-06-22 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 19. `network_dash(G, pos_dict, xx, function, color_schemes)` — Alta

- **Assinatura:** `network_dash(G, pos_dict, xx, function, color_schemes)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:137`–916 (780 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (25 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:137`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 20. `network_dash(G, pos_dict, xx, function)` — Alta

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:137`–828 (692 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (26 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:137`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 21. `network_dash(G, pos_dict, xx, function)` — Alta

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:137`–794 (658 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (25 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:137`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 22. `submitted_parameters(self, outformat)` — Alta

- **Assinatura:** `submitted_parameters(self, outformat)` (em `database`) — python, method
- **Local:** `lib/rotifer/genome/database.py:94`–103 (10 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def submitted_parameters("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 23. `available_sources(self)` — Alta

- **Assinatura:** `available_sources(self)` (em `database`) — python, method
- **Local:** `lib/rotifer/genome/database.py:105`–106 (2 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def available_sources("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 24. `fetch_one(self)` — Alta

- **Assinatura:** `fetch_one(self)` (em `clickhouse`) — python, method
- **Local:** `lib/rotifer/genome/db/clickhouse.py:154`–156 (3 linhas)
- **Última alteração (bruta):** 2019-05-23 · **excluindo ruído e importação em massa:** 2019-05-23
- **Commits que tocaram o corpo (1):** `c826e33`
- **Autores distintos (1):** kaihami
- **Introdução (pickaxe `git log -S "def fetch_one("`):** `3adf523` — 2019-03-08 — kaihami · 5 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 25. `hook_compressed_text(filename, mode, encoding)` — Alta

- **Assinatura:** `hook_compressed_text(filename, mode, encoding)` — python, function
- **Local:** `lib/rotifer/io/fileinput.py:55`–66 (12 linhas)
- **Última alteração (bruta):** 2020-05-29 · **excluindo ruído e importação em massa:** 2020-05-29
- **Commits que tocaram o corpo (1):** `e25d017`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def hook_compressed_text("`):** `7ce177f` — 2020-05-29 — Robson Francisco de Souza · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 26. `describe_modules()` — Alta

- **Assinatura:** `describe_modules()` — shell, function
- **Local:** `lib/rotifer/rbashpipe/base/00.doc.lib:6`–19 (14 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "describe_modules()"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 27. `split_array()` — Alta

- **Assinatura:** `split_array()` — shell, function
- **Local:** `lib/rotifer/rbashpipe/base/01.functions.lib:6`–8 (3 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "split_array()"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 28. `set_value_by_name()` — Alta

- **Assinatura:** `set_value_by_name()` — shell, function
- **Local:** `lib/rotifer/rbashpipe/base/01.functions.lib:43`–45 (3 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "set_value_by_name()"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 29. `make_directories()` — Alta

- **Assinatura:** `make_directories()` — shell, function
- **Local:** `lib/rotifer/rbashpipe/base/01.functions.lib:95`–103 (9 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "make_directories()"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 30. `show_options(self)` — Alta

- **Assinatura:** `show_options(self)` (em `IO`) — python, method
- **Local:** `lib/rotifer/seq/alignment.py:87`–88 (2 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `76cfb36`
- **Autores distintos (1):** kaihami
- **Introdução (pickaxe `git log -S "def show_options("`):** `76cfb36` — 2019-03-15 — kaihami · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 31. `conservation_types(self)` — Alta

- **Assinatura:** `conservation_types(self)` (em `MSA`) — python, method
- **Local:** `lib/rotifer/seq/alignment.py:638`–640 (3 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def conservation_types("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 32. `seq2colors(self)` — Alta

- **Assinatura:** `seq2colors(self)` (em `MSA`) — python, method
- **Local:** `lib/rotifer/seq/alignment.py:642`–646 (5 linhas)
- **Última alteração (bruta):** 2019-04-02 · **excluindo ruído e importação em massa:** 2019-04-02
- **Commits que tocaram o corpo (2):** `29bd8c1`, `47da1b3`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def seq2colors("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 33. `plot_logo(self, font_family, data_type, seq_type, yaxis, colorscheme, nrows, padding, draw_range, coordinate_type, draw_axis, fontfamily, debug, ax, dpi)` — Alta

- **Assinatura:** `plot_logo(self, font_family, data_type, seq_type, yaxis, colorscheme, nrows, padding, draw_range, coordinate_type, draw_axis, fontfamily, debug, ax, dpi)` (em `MSA`) — python, method
- **Local:** `lib/rotifer/seq/alignment.py:648`–674 (27 linhas)
- **Última alteração (bruta):** 2019-03-18 · **excluindo ruído e importação em massa:** 2019-03-18
- **Commits que tocaram o corpo (1):** `7d1142a`
- **Autores distintos (1):** kaihami
- **Introdução (pickaxe `git log -S "def plot_logo("`):** `7d1142a` — 2019-03-18 — kaihami · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 34. `_gap_percentage(self, matrix)` — Alta

- **Assinatura:** `_gap_percentage(self, matrix)` (em `conservation`) — python, method
- **Local:** `lib/rotifer/seq/core/core.py:83`–88 (6 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def _gap_percentage("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 35. `filter_by_size(self, col, max_length, min_length, inplace)` — Alta

- **Assinatura:** `filter_by_size(self, col, max_length, min_length, inplace)` (em `domain`) — python, method
- **Local:** `lib/rotifer/tools/search.py:438`–460 (23 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def filter_by_size("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 36. `plot_size(self, col, percentile, line_color, xlabel, ylabel, linestyle)` — Alta

- **Assinatura:** `plot_size(self, col, percentile, line_color, xlabel, ylabel, linestyle)` (em `domain`) — python, method
- **Local:** `lib/rotifer/tools/search.py:501`–523 (23 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def plot_size("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 37. `slice_sequence(self, col, start_col, end_col, inplace, col_name)` — Alta

- **Assinatura:** `slice_sequence(self, col, start_col, end_col, inplace, col_name)` (em `domain`) — python, method
- **Local:** `lib/rotifer/tools/search.py:575`–597 (23 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def slice_sequence("`):** `3adf523` — 2019-03-08 — kaihami · 7 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 38. `sub process_accotations` — Alta

- **Assinatura:** `sub process_accotations` (em `Rotifer::DBIC::AnnotationDB::Parser::genbank`) — perl, sub
- **Local:** `perl/lib/Rotifer/DBIC/AnnotationDB/Parser/genbank.pm:312`–365 (54 linhas)
- **Última alteração (bruta):** 2019-05-02 · **excluindo ruído e importação em massa:** 2019-05-02
- **Commits que tocaram o corpo (1):** `0602ce1`
- **Autores distintos (1):** Robson Francisco de Souza
- **Commits de importação em massa no intervalo (descartados):** `0602ce1`
- **Introdução (pickaxe `git log -S "sub process_accotations"`):** `0602ce1` — 2019-05-02 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo Perl instalável: scripts de laboratório fora do repo podem chamar este método diretamente — exige confirmação humana.

#### 39. `sub process_bio_searchio` — Alta

- **Assinatura:** `sub process_bio_searchio` — perl, sub
- **Local:** `bin/blast2table:129`–362 (234 linhas)
- **Última alteração (bruta):** 2025-03-11 · **excluindo ruído e importação em massa:** 2025-03-11
- **Commits que tocaram o corpo (2):** `7ee3c96`, `63adef3`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "sub process_bio_searchio"`):** `a59c36b` — 2019-11-26 — Robson Francisco de Souza · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: broken
  - bloco comentado grande (35 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 40. `hierarchy_to_dataframe(x)` — Alta

- **Assinatura:** `hierarchy_to_dataframe(x)` — python, function
- **Local:** `lib/rotifer/core/functions.py:385`–397 (13 linhas)
- **Última alteração (bruta):** 2022-06-10 · **excluindo ruído e importação em massa:** 2022-06-10
- **Commits que tocaram o corpo (1):** `32f427f`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def hierarchy_to_dataframe("`):** `32f427f` — 2022-06-09 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 41. `aln_fig_style(i, polish_aln, consensus, annotations, remove_gaps, adjust_coordinates, font_size)` — Alta

- **Assinatura:** `aln_fig_style(i, polish_aln, consensus, annotations, remove_gaps, adjust_coordinates, font_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:4192`–4326 (135 linhas)
- **Última alteração (bruta):** 2025-05-12 · **excluindo ruído e importação em massa:** 2025-05-12
- **Commits que tocaram o corpo (1):** `bfc33a9`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def aln_fig_style("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: todo
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 42. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Alta

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:12`–110 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:12`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 43. `update_stylesheet(font_size, outline_size, edge_size)` — Alta

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:391`–418 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:389`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 44. `update_input_boxes(node_data, toggle)` — Alta

- **Assinatura:** `update_input_boxes(node_data, toggle)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:425`–432 (8 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_input_boxes("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:465`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 45. `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, selected_scheme, elements, last_click, color_cycle, last_node_click, node_color_cycle, color_schemes, pos_dict_data, serial_graph)` — Alta

- **Assinatura:** `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, selected_scheme, elements, last_click, color_cycle, last_node_click, node_color_cycle, color_schemes, pos_dict_data, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:467`–634 (168 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_elements_or_cycle("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 7 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:507`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 46. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` — Alta

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:650`–664 (15 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:676`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 47. `print_stored_data(n_clicks, pos_dict)` — Alta

- **Assinatura:** `print_stored_data(n_clicks, pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:675`–687 (13 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:695`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 48. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Alta

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:701`–710 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:710`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 49. `update_group_dropdown(color_schemes_dicts, colorscheme)` — Alta

- **Assinatura:** `update_group_dropdown(color_schemes_dicts, colorscheme)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:719`–723 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_group_dropdown("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 5 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions2.py:713 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 50. `add_group(n_clicks, new_key, new_color, color_schemes_dicts, colorscheme)` — Alta

- **Assinatura:** `add_group(n_clicks, new_key, new_color, color_schemes_dicts, colorscheme)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:734`–736 (3 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def add_group("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:752`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 51. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph, function_data)` — Alta

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph, function_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:752`–767 (16 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:771`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 52. `modify_edge_and_store_current(add_clicks, remove_clicks, source, target, graph_data, trigger_count, elements, pos_dict)` — Alta

- **Assinatura:** `modify_edge_and_store_current(add_clicks, remove_clicks, source, target, graph_data, trigger_count, elements, pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:786`–815 (30 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def modify_edge_and_store_current("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório — todas na mesma família de arquivo (original + cópias de backup), portanto a ambiguidade é evidência de duplicação, não de uso
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 53. `populate_colorscheme_dropdown(color_data)` — Alta

- **Assinatura:** `populate_colorscheme_dropdown(color_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:822`–826 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def populate_colorscheme_dropdown("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório — todas na mesma família de arquivo (original + cópias de backup), portanto a ambiguidade é evidência de duplicação, não de uso
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:821`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 54. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Alta

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:953`–980 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:870`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 55. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Alta

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:982`–1004 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:899`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 56. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Alta

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:1006`–1094 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:923`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 57. `update_group_dropdown(xx_data)` — Alta

- **Assinatura:** `update_group_dropdown(xx_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:713`–715 (3 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_group_dropdown("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 5 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions2.py:713 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 58. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Alta

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:12`–110 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:12`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 59. `update_stylesheet(font_size, outline_size, edge_size)` — Alta

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:387`–414 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:389`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 60. `update_input_boxes(node_data, toggle)` — Alta

- **Assinatura:** `update_input_boxes(node_data, toggle)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:463`–470 (8 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_input_boxes("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:465`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 61. `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` — Alta

- **Assinatura:** `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:505`–659 (155 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_elements_or_cycle("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 7 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:507`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 62. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` — Alta

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:675`–685 (11 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:676`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 63. `print_stored_data(n_clicks, toprint)` — Alta

- **Assinatura:** `print_stored_data(n_clicks, toprint)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:694`–697 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:695`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 64. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Alta

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:709`–718 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:710`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 65. `update_group_dropdown(xx_data, function_data)` — Alta

- **Assinatura:** `update_group_dropdown(xx_data, function_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:727`–730 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_group_dropdown("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 5 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions2.py:713 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 66. `add_group(n_clicks, new_key, new_color, current_xx)` — Alta

- **Assinatura:** `add_group(n_clicks, new_key, new_color, current_xx)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:740`–744 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def add_group("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:752`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 67. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph, saved_func_dict)` — Alta

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph, saved_func_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:760`–774 (15 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:771`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 68. `modify_edge_and_store_current(add_clicks, remove_clicks, source, target, graph_data, trigger_count, elements, pos_dict)` — Alta

- **Assinatura:** `modify_edge_and_store_current(add_clicks, remove_clicks, source, target, graph_data, trigger_count, elements, pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:793`–822 (30 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def modify_edge_and_store_current("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório — todas na mesma família de arquivo (original + cópias de backup), portanto a ambiguidade é evidência de duplicação, não de uso
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 69. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Alta

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:865`–892 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:870`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 70. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Alta

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:894`–916 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:899`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 71. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Alta

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:918`–1006 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:923`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 72. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Alta

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:12`–110 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:12`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 73. `update_stylesheet(font_size, outline_size, edge_size)` — Alta

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:375`–402 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:389`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 74. `update_input_boxes(node_data, toggle)` — Alta

- **Assinatura:** `update_input_boxes(node_data, toggle)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:451`–458 (8 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_input_boxes("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:465`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 75. `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` — Alta

- **Assinatura:** `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:493`–646 (154 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_elements_or_cycle("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 7 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:507`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 76. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` — Alta

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:662`–672 (11 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:676`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 77. `print_stored_data(n_clicks, toprint)` — Alta

- **Assinatura:** `print_stored_data(n_clicks, toprint)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:681`–684 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:695`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 78. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Alta

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:696`–705 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:710`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 79. `update_group_dropdown(xx_data)` — Alta

- **Assinatura:** `update_group_dropdown(xx_data)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:713`–715 (3 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_group_dropdown("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 5 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions2.py:713 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 80. `add_group(n_clicks, new_key, new_color, current_xx)` — Alta

- **Assinatura:** `add_group(n_clicks, new_key, new_color, current_xx)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:724`–728 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def add_group("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:752`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 81. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph)` — Alta

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:743`–757 (15 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:771`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 82. `modify_edge(add_clicks, remove_clicks, source, target, graph_data, trigger_count)` — Alta

- **Assinatura:** `modify_edge(add_clicks, remove_clicks, source, target, graph_data, trigger_count)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:771`–786 (16 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def modify_edge("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 4 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:799`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 83. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Alta

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:831`–858 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:870`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 84. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Alta

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:860`–882 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:899`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 85. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Alta

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:884`–972 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:923`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 86. `optimize_memory_usage(df)` — Alta

- **Assinatura:** `optimize_memory_usage(df)` — python, function
- **Local:** `lib/rotifer/pandas/functions.py:97`–140 (44 linhas)
- **Última alteração (bruta):** 2021-12-21 · **excluindo ruído e importação em massa:** 2021-12-21
- **Commits que tocaram o corpo (1):** `0efbcc7`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def optimize_memory_usage("`):** `0efbcc7` — 2021-12-21 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 87. `print_everything(max_rows, max_columns, max_colwidth, width, verbose)` — Alta

- **Assinatura:** `print_everything(max_rows, max_columns, max_colwidth, width, verbose)` — python, function
- **Local:** `lib/rotifer/pandas/functions.py:142`–170 (29 linhas)
- **Última alteração (bruta):** 2021-12-21 · **excluindo ruído e importação em massa:** 2021-12-21
- **Commits que tocaram o corpo (2):** `0efbcc7`, `62d529e`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def print_everything("`):** `62d529e` — 2020-07-30 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 88. `show_display_options()` — Alta

- **Assinatura:** `show_display_options()` — python, function
- **Local:** `lib/rotifer/pandas/functions.py:172`–177 (6 linhas)
- **Última alteração (bruta):** 2021-12-21 · **excluindo ruído e importação em massa:** 2021-12-21
- **Commits que tocaram o corpo (1):** `0efbcc7`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def show_display_options("`):** `0efbcc7` — 2021-12-21 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 89. `calculate_positions(G, iterations)` — Alta

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:112`–121 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:112`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 90. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Alta

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:123`–134 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:123`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 91. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Alta

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:919`–951 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:836`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 92. `calculate_positions(G, iterations)` — Alta

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:112`–121 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:112`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 93. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Alta

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:123`–134 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:123`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 94. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Alta

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.quebrada.py:831`–863 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:836`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 95. `calculate_positions(G, iterations)` — Alta

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:112`–121 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:112`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 96. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Alta

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:123`–134 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:123`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 97. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Alta

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.working_20250627.py:797`–829 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:836`
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

---

## 2. Média confiança — depreciar com aviso antes de remover — 72 funções

Dois eixos convergem, o terceiro é fraco ou ausente. O padrão recomendado é emitir `DeprecationWarning` (ou comentário equivalente em Perl/Shell) por um ciclo antes da remoção.

#### 1. `arch_2_svg(df, length, font_size)` — Média

- **Assinatura:** `arch_2_svg(df, length, font_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/draw.py:1`–127 (127 linhas)
- **Última alteração (bruta):** 2024-06-09 · **excluindo ruído e importação em massa:** 2024-06-09
- **Commits que tocaram o corpo (4):** `2f40ea7`, `e5165d6`, `f67cc8c`, `6cd098c`
- **Autores distintos (2):** Robson Francisco de Souza, rodolfoar
- **Introdução (pickaxe `git log -S "def arch_2_svg("`):** `e5165d6` — 2024-03-06 — rodolfoar
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (24 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 2. `uniref50_to_ncbi(uniref50_strings)` — Média

- **Assinatura:** `uniref50_to_ncbi(uniref50_strings)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1066`–1087 (22 linhas)
- **Última alteração (bruta):** 2024-06-21 · **excluindo ruído e importação em massa:** 2024-06-21
- **Commits que tocaram o corpo (1):** `dbc2d64`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def uniref50_to_ncbi("`):** `dbc2d64` — 2024-06-21 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 3. `uniref50_to_clusters(uniref50_strings)` — Média

- **Assinatura:** `uniref50_to_clusters(uniref50_strings)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1089`–1123 (35 linhas)
- **Última alteração (bruta):** 2024-06-21 · **excluindo ruído e importação em massa:** 2024-06-21
- **Commits que tocaram o corpo (1):** `dbc2d64`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def uniref50_to_clusters("`):** `dbc2d64` — 2024-06-21 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 4. `pid2uniref50(pids)` — Média

- **Assinatura:** `pid2uniref50(pids)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1125`–1152 (28 linhas)
- **Última alteração (bruta):** 2024-06-21 · **excluindo ruído e importação em massa:** 2024-06-21
- **Commits que tocaram o corpo (1):** `dbc2d64`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def pid2uniref50("`):** `dbc2d64` — 2024-06-21 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 5. `operon_fig_bkp(df, domain_rename, output_file, domain_column, query_same_direction, color_dict, top_domains, sort_by, height, f, fontsize, query_asterix, query_color, check_duplicates, light_palette, reverse_annotaion)` — Média

- **Assinatura:** `operon_fig_bkp(df, domain_rename, output_file, domain_column, query_same_direction, color_dict, top_domains, sort_by, height, f, fontsize, query_asterix, query_color, check_duplicates, light_palette, reverse_annotaion)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1885`–2077 (193 linhas)
- **Última alteração (bruta):** 2024-12-05 · **excluindo ruído e importação em massa:** 2024-12-05
- **Commits que tocaram o corpo (12):** `672e490`, `4329acd`, `82d3061`, `56e9a81`, `c9d3bba`, `7859d49`, `0d76404`, `db67e8a`…
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def operon_fig_bkp("`):** `56e9a81` — 2024-10-24 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (13 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 6. `clean_dict(d)` — Média

- **Assinatura:** `clean_dict(d)` (em `network_dash.print_stored_data`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:676`–683 (8 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clean_dict("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:685`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 7. `update_graph(selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value)` — Média

- **Assinatura:** `update_graph(selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:349`–377 (29 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_graph("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 6** — `lib/rotifer/devel/alpha/net_functions.bkp2:356`, `lib/rotifer/devel/alpha/net_functions.good.bkp:414`, `lib/rotifer/devel/alpha/net_functions.py:428`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:426`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:414`, `lib/rotifer/devel/alpha/net_functions2.py:414`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 8. `network_dash(G, pos_dict, xx, function)` — Média

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:136`–638 (503 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (16 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 9. `network_dash(G, pos_dict, xx, function)` — Média

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:137`–794 (658 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (25 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 10. `get_some_info()` — Média

- **Assinatura:** `get_some_info()` — python, function
- **Local:** `lib/rotifer/core/functions.py:29`–76 (48 linhas)
- **Última alteração (bruta):** 2024-02-10 · **excluindo ruído e importação em massa:** 2024-02-10
- **Commits que tocaram o corpo (1):** `05bbe11`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def get_some_info("`):** `05bbe11` — 2024-02-10 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (12 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 11. `SequenceCollection` — Média

- **Assinatura:** `SequenceCollection` — python, class
- **Local:** `lib/rotifer/devel/alpha/collection.py:25`–138 (114 linhas)
- **Última alteração (bruta):** 2023-10-11 · **excluindo ruído e importação em massa:** 2023-10-11
- **Commits que tocaram o corpo (5):** `45d8511`, `ed2ce85`, `a0d0961`, `0b2844f`, `f12f77b`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "class SequenceCollection("`):** `ed2ce85` — 2023-06-01 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 12. `get_correspondent_position(seqobj_source, seqobj_target, pid, position)` — Média

- **Assinatura:** `get_correspondent_position(seqobj_source, seqobj_target, pid, position)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1001`–1018 (18 linhas)
- **Última alteração (bruta):** 2024-05-17 · **excluindo ruído e importação em massa:** 2024-05-17
- **Commits que tocaram o corpo (3):** `f9c4957`, `ce56353`, `fdd32b1`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_correspondent_position("`):** `f9c4957` — 2024-04-09 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 13. `read_predicted_topologies(file)` — Média

- **Assinatura:** `read_predicted_topologies(file)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1021`–1045 (25 linhas)
- **Última alteração (bruta):** 2024-05-17 · **excluindo ruído e importação em massa:** 2024-05-17
- **Commits que tocaram o corpo (1):** `ce56353`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def read_predicted_topologies("`):** `ce56353` — 2024-05-17 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 14. `uniref50_add_info(seqobj)` — Média

- **Assinatura:** `uniref50_add_info(seqobj)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1047`–1063 (17 linhas)
- **Última alteração (bruta):** 2024-05-28 · **excluindo ruído e importação em massa:** 2024-05-28
- **Commits que tocaram o corpo (2):** `f64855b`, `da74d29`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def uniref50_add_info("`):** `f64855b` — 2024-05-24 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 15. `mview(seqobj, output, background, consensus, organism, find)` — Média

- **Assinatura:** `mview(seqobj, output, background, consensus, organism, find)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1157`–1192 (36 linhas)
- **Última alteração (bruta):** 2024-06-21 · **excluindo ruído e importação em massa:** 2024-06-21
- **Commits que tocaram o corpo (1):** `dbc2d64`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def mview("`):** `dbc2d64` — 2024-06-21 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 16. `uniprot_to_ncbi(uniprot_id)` — Média

- **Assinatura:** `uniprot_to_ncbi(uniprot_id)` — python, function
- **Local:** `lib/rotifer/devel/alpha/idmapping_uniprot.py:222`–237 (16 linhas)
- **Última alteração (bruta):** 2024-05-29 · **excluindo ruído e importação em massa:** 2024-05-29
- **Commits que tocaram o corpo (1):** `3c34b64`
- **Autores distintos (1):** Rodolfo Alvarenga Ribeiro
- **Introdução (pickaxe `git log -S "def uniprot_to_ncbi("`):** `3c34b64` — 2024-05-29 — Rodolfo Alvarenga Ribeiro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - url externa: https://rest.uniprot.org/uniprotkb/{uniprot_id}.json
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 17. `filter_neighbors_plus(ndf, pids, reqdom, customdoms, after, before, max_distance, genome_protein_fasta, models_path, annotate, mode, seed, patience, max_extend)` — Média

- **Assinatura:** `filter_neighbors_plus(ndf, pids, reqdom, customdoms, after, before, max_distance, genome_protein_fasta, models_path, annotate, mode, seed, patience, max_extend)` — python, function
- **Local:** `lib/rotifer/devel/alpha/malu.py:499`–657 (159 linhas)
- **Última alteração (bruta):** 2026-09-07 · **excluindo ruído e importação em massa:** 2026-09-07
- **Commits que tocaram o corpo (2):** `2b05777`, `7b9badf`
- **Autores distintos (1):** Maria Luiza Andreani
- **Introdução (pickaxe `git log -S "def filter_neighbors_plus("`):** `2b05777` — 2026-08-13 — Maria Luiza Andreani
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/epsoares.py:1340`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: broken
  - bloco comentado grande (13 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Maria Luiza Andreani`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 18. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Média

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:11`–109 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 19. `network_dash(G, pos_dict, xx, function)` — Média

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:136`–461 (326 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 20. `update_stylesheet(font_size, outline_size, edge_size)` — Média

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:310`–337 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 21. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements)` — Média

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:390`–400 (11 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 22. `print_stored_data(n_clicks, toprint)` — Média

- **Assinatura:** `print_stored_data(n_clicks, toprint)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:409`–412 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 23. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Média

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:424`–433 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 24. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections)` — Média

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:443`–456 (14 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 25. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Média

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:498`–525 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 26. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Média

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:527`–549 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 27. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Média

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:551`–639 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 28. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Média

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:11`–109 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 29. `update_stylesheet(font_size, outline_size, edge_size)` — Média

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:317`–344 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 30. `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, elements, last_click, color_cycle, last_node_click, node_color_cycle)` — Média

- **Assinatura:** `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, elements, last_click, color_cycle, last_node_click, node_color_cycle)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:413`–549 (137 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_elements_or_cycle("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 7 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 31. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements)` — Média

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:562`–572 (11 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 32. `print_stored_data(n_clicks, toprint)` — Média

- **Assinatura:** `print_stored_data(n_clicks, toprint)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:581`–584 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 33. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Média

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:596`–605 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 34. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections)` — Média

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:615`–628 (14 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 35. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Média

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:675`–702 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 36. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Média

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:704`–726 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 37. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Média

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:728`–816 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 38. `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — Média

- **Assinatura:** `network_annotation(networkdf, node_color, node_shape, second_shape, second_shape_list, node_outline, node_outline_color, node_outline_list, edge_color, unidirected_graph)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:12`–110 (99 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_annotation("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 39. `update_stylesheet(font_size, outline_size, edge_size)` — Média

- **Assinatura:** `update_stylesheet(font_size, outline_size, edge_size)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:375`–402 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_stylesheet("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 40. `update_input_boxes(node_data, toggle)` — Média

- **Assinatura:** `update_input_boxes(node_data, toggle)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:451`–458 (8 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_input_boxes("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 41. `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` — Média

- **Assinatura:** `update_elements_or_cycle(edge_data, node_data, selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value, pos_dict_input, dummy_call, elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:493`–646 (154 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_elements_or_cycle("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 7 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 42. `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` — Média

- **Assinatura:** `save_positions(n_clicks, text, x_scale, y_scale, scale, n_elements, current_pos_dict)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:662`–672 (11 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def save_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 43. `print_stored_data(n_clicks, toprint)` — Média

- **Assinatura:** `print_stored_data(n_clicks, toprint)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:681`–684 (4 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def print_stored_data("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 44. `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` — Média

- **Assinatura:** `get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:696`–705 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_image("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 45. `add_group(n_clicks, new_key, new_color, current_xx)` — Média

- **Assinatura:** `add_group(n_clicks, new_key, new_color, current_xx)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:724`–728 (5 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def add_group("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 6 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 46. `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph)` — Média

- **Assinatura:** `update_layouts(n_clicks, groups, Spring_lay_iteractions, min_connections, pos_dict, serial_graph)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:743`–757 (15 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def update_layouts("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 47. `modify_edge(add_clicks, remove_clicks, source, target, graph_data, trigger_count)` — Média

- **Assinatura:** `modify_edge(add_clicks, remove_clicks, source, target, graph_data, trigger_count)` (em `network_dash`) — python, method
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:771`–786 (16 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def modify_edge("`):** `1739875` — 2025-08-01 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 4 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 48. `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — Média

- **Assinatura:** `clique_jaccard_community(G, min_jaccard_index, resolution_parameter, min_clique_size)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:831`–858 (28 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def clique_jaccard_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 49. `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — Média

- **Assinatura:** `subgraph_stats(net, subnet, theme_dict, theme_to_stats, alternative, return_df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:860`–882 (23 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def subgraph_stats("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 50. `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — Média

- **Assinatura:** `get_relations(domdf, dom, min_connections, self_relation, to_yaml, iteration, filter_rename_yaml)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:884`–972 (89 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_relations("`):** `e21006c` — 2024-11-07 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 9 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.

#### 51. `check_rneighbors(acc_file, rneighbors_file)` — Média

- **Assinatura:** `check_rneighbors(acc_file, rneighbors_file)` (em `rsearch`) — python, method
- **Local:** `bin/rsearch:251`–252 (2 linhas)
- **Última alteração (bruta):** 2019-03-08 · **excluindo ruído e importação em massa:** 2019-03-08
- **Commits que tocaram o corpo (1):** `3adf523`
- **Autores distintos (1):** kaihami
- **Commits de importação em massa no intervalo (descartados):** `3adf523`
- **Introdução (pickaxe `git log -S "def check_rneighbors("`):** `3adf523` — 2019-03-08 — kaihami · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 52. `check_rneighbors(acc_file, rneighbors_file)` — Média

- **Assinatura:** `check_rneighbors(acc_file, rneighbors_file)` — python, function
- **Local:** `bin/rsearch2:226`–227 (2 linhas)
- **Última alteração (bruta):** 2019-03-08 · **excluindo ruído e importação em massa:** 2019-03-08
- **Commits que tocaram o corpo (1):** `3adf523`
- **Autores distintos (1):** kaihami
- **Commits de importação em massa no intervalo (descartados):** `3adf523`
- **Introdução (pickaxe `git log -S "def check_rneighbors("`):** `3adf523` — 2019-03-08 — kaihami · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 53. `check_rneighbors(acc_file, rneighbors_file)` — Média

- **Assinatura:** `check_rneighbors(acc_file, rneighbors_file)` — python, function
- **Local:** `bin/rsearch_old:243`–244 (2 linhas)
- **Última alteração (bruta):** 2019-03-08 · **excluindo ruído e importação em massa:** 2019-03-08
- **Commits que tocaram o corpo (1):** `3adf523`
- **Autores distintos (1):** kaihami
- **Commits de importação em massa no intervalo (descartados):** `3adf523`
- **Introdução (pickaxe `git log -S "def check_rneighbors("`):** `3adf523` — 2019-03-08 — kaihami · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `bin/rsearch:251`
- **Risco de quebra retroativa:** Confinado ao próprio script: a função não é alcançável de fora do arquivo, então o risco é apenas de algum caminho de execução interno não exercitado.

#### 54. `dict_options_dev` — Média

- **Assinatura:** `dict_options_dev` (em `action`) — python, class
- **Local:** `lib/rotifer/core/cli.py:753`–771 (19 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "class dict_options_dev("`):** `3adf523` — 2019-03-08 — kaihami · 5 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/core/dev/cli.py:707 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 55. `dict_options_dev` — Média

- **Assinatura:** `dict_options_dev` (em `action`) — python, class
- **Local:** `lib/rotifer/core/dev/cli.py:707`–725 (19 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "class dict_options_dev("`):** `3adf523` — 2019-03-08 — kaihami · 5 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/core/cli.py:753 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 56. `CoreMethods` — Média

- **Assinatura:** `CoreMethods` — python, class
- **Local:** `lib/rotifer/core/methods.py:4`–6 (3 linhas)
- **Última alteração (bruta):** 2024-02-26 · **excluindo ruído e importação em massa:** 2024-02-26
- **Commits que tocaram o corpo (1):** `f8b3743`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe):** [sem evidência]
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 57. `FileCollection` — Média

- **Assinatura:** `FileCollection` — python, class
- **Local:** `lib/rotifer/db/local/core.py:18`–90 (73 linhas)
- **Última alteração (bruta):** 2023-12-14 · **excluindo ruído e importação em massa:** 2023-12-14
- **Commits que tocaram o corpo (1):** `db6c748`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "class FileCollection("`):** `db6c748` — 2023-12-14 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 58. `update_database(self)` — Média

- **Assinatura:** `update_database(self)` (em `TaxonomyCursor`) — python, method
- **Local:** `lib/rotifer/db/local/ete3.py:33`–42 (10 linhas)
- **Última alteração (bruta):** 2022-12-20 · **excluindo ruído e importação em massa:** 2022-12-20
- **Commits que tocaram o corpo (2):** `50da12e`, `49d6c04`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def update_database("`):** `50da12e` — 2022-11-16 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/trsantos.py:302`, `lib/rotifer/devel/alpha/trsantos.py:475`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 59. `getids2(self, obj, *args, **kwargs)` — Média

- **Assinatura:** `getids2(self, obj, *args, **kwargs)` (em `GeneNeighborhoodCursor`) — python, method
- **Local:** `lib/rotifer/db/ncbi/entrez.py:485`–493 (9 linhas)
- **Última alteração (bruta):** 2023-01-23 · **excluindo ruído e importação em massa:** 2023-01-23
- **Commits que tocaram o corpo (2):** `db89521`, `b00dcd2`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def getids2("`):** `b00dcd2` — 2023-01-23 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 60. `cluster2aln(group_cluster, df, esl_index_file, grouper, redundancy_cluster, align_method, query, cpu)` — Média

- **Assinatura:** `cluster2aln(group_cluster, df, esl_index_file, grouper, redundancy_cluster, align_method, query, cpu)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:65`–115 (51 linhas)
- **Última alteração (bruta):** 2022-10-18 · **excluindo ruído e importação em massa:** 2022-10-18
- **Commits que tocaram o corpo (6):** `db472c1`, `2591006`, `c42772d`, `d3614dc`, `86cc0b7`, `4c7073c`
- **Autores distintos (3):** Gianlucca Gonçalves Nicastro, Robson Francisco de Souza, rodolfoar
- **Introdução (pickaxe `git log -S "def cluster2aln("`):** `2591006` — 2022-05-25 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 61. `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — Média

- **Assinatura:** `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:124`–193 (70 linhas)
- **Última alteração (bruta):** 2024-06-09 · **excluindo ruído e importação em massa:** 2024-06-09
- **Commits que tocaram o corpo (5):** `db472c1`, `a4257e2`, `0c86ab9`, `d24dddc`, `6cd098c`
- **Autores distintos (2):** Gianlucca Gonçalves Nicastro, rodolfoar
- **Introdução (pickaxe `git log -S "def cluster_Co_occurrence("`):** `0c86ab9` — 2022-04-14 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: xxx,xxxx
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 62. `hmmsearch_full2pandas(file, error_lines, keep_threshold)` — Média

- **Assinatura:** `hmmsearch_full2pandas(file, error_lines, keep_threshold)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:338`–385 (48 linhas)
- **Última alteração (bruta):** 2022-10-18 · **excluindo ruído e importação em massa:** 2022-10-18
- **Commits que tocaram o corpo (4):** `db472c1`, `6f5b03c`, `6d86092`, `a80a16d`
- **Autores distintos (3):** Gianlucca G. Nicastro, Gianlucca Gonçalves Nicastro, Gianlucca Nicastro
- **Introdução (pickaxe `git log -S "def hmmsearch_full2pandas("`):** `db472c1` — 2022-10-17 — Gianlucca Gonçalves Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 63. `hhr_to_aln(seqobj, hhr, database)` — Média

- **Assinatura:** `hhr_to_aln(seqobj, hhr, database)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:388`–481 (94 linhas)
- **Última alteração (bruta):** 2022-10-18 · **excluindo ruído e importação em massa:** 2022-10-18
- **Commits que tocaram o corpo (4):** `df6d548`, `db472c1`, `f88f1d5`, `aedec01`
- **Autores distintos (2):** Gianlucca G. Nicastro, Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def hhr_to_aln("`):** `df6d548` — 2022-08-02 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 64. `add_arch_to_seqobj(seqobj, db, cpu)` — Média

- **Assinatura:** `add_arch_to_seqobj(seqobj, db, cpu)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:771`–796 (26 linhas)
- **Última alteração (bruta):** 2022-10-14 · **excluindo ruído e importação em massa:** 2022-10-14
- **Commits que tocaram o corpo (2):** `eae51c2`, `8124a4b`
- **Autores distintos (1):** Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def add_arch_to_seqobj("`):** `eae51c2` — 2022-10-14 — Gianlucca Gonçalves Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 65. `scaled_repeat_region_svg(x_start, x_end, lane_top, lane_bottom, connector_y, repeat_start, repeat_end, repeat_strand, flip, clamp_min, clamp_max)` — Média

- **Assinatura:** `scaled_repeat_region_svg(x_start, x_end, lane_top, lane_bottom, connector_y, repeat_start, repeat_end, repeat_strand, flip, clamp_min, clamp_max)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:1381`–1446 (66 linhas)
- **Última alteração (bruta):** 2026-09-09 · **excluindo ruído e importação em massa:** 2026-09-09
- **Commits que tocaram o corpo (4):** `47d9d71`, `dc24cee`, `2bd7e3f`, `7a525ef`
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def scaled_repeat_region_svg("`):** `2bd7e3f` — 2026-08-28 — EdwardNSeven
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 3** — `lib/rotifer/devel/alpha/igem.py:101`, `lib/rotifer/devel/alpha/igem.py:2954`, `lib/rotifer/devel/alpha/igem.py:3171`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: no longer
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 66. `annotate_community(df, column, reference)` — Média

- **Assinatura:** `annotate_community(df, column, reference)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:17`–33 (17 linhas)
- **Última alteração (bruta):** 2022-08-08 · **excluindo ruído e importação em massa:** 2022-08-08
- **Commits que tocaram o corpo (1):** `ff827df`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def annotate_community("`):** `ff827df` — 2022-08-08 — rodolfoar · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:17 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 67. `dataframe_to_community(df, source, target, weight, directed)` — Média

- **Assinatura:** `dataframe_to_community(df, source, target, weight, directed)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:35`–61 (27 linhas)
- **Última alteração (bruta):** 2022-08-08 · **excluindo ruído e importação em massa:** 2022-08-08
- **Commits que tocaram o corpo (3):** `7a35537`, `d3df041`, `ff827df`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def dataframe_to_community("`):** `7a35537` — 2022-08-08 — rodolfoar · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:35 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 68. `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — Média

- **Assinatura:** `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:272`–341 (70 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def cluster_Co_occurrence("`):** `0c86ab9` — 2022-04-14 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: xxx,xxxx
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.

#### 69. `parse_hhr(indir, suffix)` — Média

- **Assinatura:** `parse_hhr(indir, suffix)` — python, function
- **Local:** `lib/rotifer/io/hhsuite.py:46`–74 (29 linhas)
- **Última alteração (bruta):** 2024-03-16 · **excluindo ruído e importação em massa:** 2024-03-16
- **Commits que tocaram o corpo (3):** `8cefc6d`, `5e64a0c`, `2582bdb`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def parse_hhr("`):** `8cefc6d` — 2022-01-05 — Robson Francisco de Souza · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 70. `seq_len(self, col)` — Média

- **Assinatura:** `seq_len(self, col)` (em `protein`) — python, method
- **Local:** `lib/rotifer/seq/seq.py:57`–67 (11 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def seq_len("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/tools/search.py:599 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 71. `seq_len(self, col, col_name, inplace)` — Média

- **Assinatura:** `seq_len(self, col, col_name, inplace)` (em `domain`) — python, method
- **Local:** `lib/rotifer/tools/search.py:599`–613 (15 linhas)
- **Última alteração (bruta):** 2019-03-15 · **excluindo ruído e importação em massa:** 2019-03-15
- **Commits que tocaram o corpo (1):** `29bd8c1`
- **Autores distintos (1):** kaihami
- **Commits de ruído no intervalo (descartados):** `29bd8c1`
- **Introdução (pickaxe `git log -S "def seq_len("`):** `3adf523` — 2019-03-08 — kaihami · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/seq/seq.py:57 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

#### 72. `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — Média

- **Assinatura:** `cluster_Co_occurrence(df, count, freq_cutoff, only_query, annotation)` — python, function
- **Local:** `lib/rotifer/devel/beta/blast.py:272`–341 (70 linhas)
- **Última alteração (bruta):** 2025-07-14 · **excluindo ruído e importação em massa:** 2025-07-14
- **Commits que tocaram o corpo (1):** `a1d6863`
- **Autores distintos (1):** Eduardo Pereira Soares
- **Introdução (pickaxe `git log -S "def cluster_Co_occurrence("`):** `0c86ab9` — 2022-04-14 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: xxx,xxxx
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.

---

## 3. Baixa confiança — investigar — 61 funções

Evidência conflitante. Cada ficha termina com a pergunta exata que precisa ser respondida e por quem.

#### 1. `igem_pipeline(genome_annotation, genome_format, genome_protein_fasta, genome_nucleotide_fasta, models_path, search_models, hmmsearch_score_filter, hmmsearch_evalue_filter, return_hmmscan, after, before, run_fimo, meme_file, return_fimo, make_figure, output_report, repeat_max_distance, repeat_min_spacing, repeat_max_spacing, min_repeats, color_dict, domain_dict, seed, patience, max_distance, max_extend, domains_filter, organism, add_sequences, normalize_orientation, filter_columns)` — Baixa

- **Assinatura:** `igem_pipeline(genome_annotation, genome_format, genome_protein_fasta, genome_nucleotide_fasta, models_path, search_models, hmmsearch_score_filter, hmmsearch_evalue_filter, return_hmmscan, after, before, run_fimo, meme_file, return_fimo, make_figure, output_report, repeat_max_distance, repeat_min_spacing, repeat_max_spacing, min_repeats, color_dict, domain_dict, seed, patience, max_distance, max_extend, domains_filter, organism, add_sequences, normalize_orientation, filter_columns)` — python, function
- **Local:** `lib/rotifer/devel/alpha/epsoares.py:1311`–1400 (90 linhas)
- **Última alteração (bruta):** 2026-09-10 · **excluindo ruído e importação em massa:** 2026-09-10
- **Commits que tocaram o corpo (20):** `2cc263d`, `dc24cee`, `de77e4e`, `0e8d85f`, `0d173a6`, `fe4a5d4`, `27a57b0`, `2bd7e3f`…
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def igem_pipeline("`):** `de77e4e` — 2026-06-16 — Eduardo Pereira Soares
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /databases/pfam/Pfam-A.hmm; /home/leep/epsoares/projects/igem/2026/data/all_models.hmm; /home/leep/epsoares/projects/igem/2026/data/heptarepeats2.meme
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Eduardo Pereira Soares`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Eduardo Pereira Soares

#### 2. `draw_architecture(input_file, file, id_list, size_median, not_spaces, domain_rename, db, evalue, color_dict, shape_dict, round_dict)` — Baixa

- **Assinatura:** `draw_architecture(input_file, file, id_list, size_median, not_spaces, domain_rename, db, evalue, color_dict, shape_dict, round_dict)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1382`–1646 (265 linhas)
- **Última alteração (bruta):** 2026-02-26 · **excluindo ruído e importação em massa:** 2026-02-26
- **Commits que tocaram o corpo (8):** `90dddc5`, `a93375d`, `ede3e25`, `b8a044a`, `c9d3bba`, `558f52c`, `0d76404`, `82d3061`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def draw_architecture("`):** `90dddc5` — 2024-09-09 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (33 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Gianlucca G. Nicastro

#### 3. `mmseqs_easy_search(acc, cpu, aln, max_out, sensitivity, iterations, db, columns)` — Baixa

- **Assinatura:** `mmseqs_easy_search(acc, cpu, aln, max_out, sensitivity, iterations, db, columns)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:3982`–4056 (75 linhas)
- **Última alteração (bruta):** 2025-05-21 · **excluindo ruído e importação em massa:** 2025-05-21
- **Commits que tocaram o corpo (2):** `bfc33a9`, `8495bbd`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def mmseqs_easy_search("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/data/fadb/tmp/nr.50.mmseqs.db
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Gianlucca G. Nicastro

#### 4. `run_dommain_annotation_farm(input_file, cores, batch_size, profiledb, pfam, hmm_db_profiledb, hmm_db_pfam, jobs)` — Baixa

- **Assinatura:** `run_dommain_annotation_farm(input_file, cores, batch_size, profiledb, pfam, hmm_db_profiledb, hmm_db_pfam, jobs)` — python, function
- **Local:** `lib/rotifer/devel/alpha/snakemake.py:270`–300 (31 linhas)
- **Última alteração (bruta):** 2025-05-23 · **excluindo ruído e importação em massa:** 2025-05-23
- **Commits que tocaram o corpo (2):** `8619346`, `8495bbd`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def run_dommain_annotation_farm("`):** `8619346` — 2025-05-23 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/data/rpsdb/allprofiles; /netmnt/vast01/cbb01/proteinworld/data/rpsdb/pwld_new_pfam
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Gianlucca G. Nicastro

#### 5. `run_dommain_annotation_farm_from_acc(input_file, cores, batch_size, profiledb, pfam, hmm_db_profiledb, hmm_db_pfam, delete_tmpdir, existing_tmpdir, jobs)` — Baixa

- **Assinatura:** `run_dommain_annotation_farm_from_acc(input_file, cores, batch_size, profiledb, pfam, hmm_db_profiledb, hmm_db_pfam, delete_tmpdir, existing_tmpdir, jobs)` — python, function
- **Local:** `lib/rotifer/devel/alpha/snakemake.py:306`–340 (35 linhas)
- **Última alteração (bruta):** 2026-03-02 · **excluindo ruído e importação em massa:** 2026-03-02
- **Commits que tocaram o corpo (2):** `705301e`, `42a897b`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def run_dommain_annotation_farm_from_acc("`):** `705301e` — 2025-06-13 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/data/rpsdb/allprofiles; /netmnt/vast01/cbb01/proteinworld/data/rpsdb/pwld_new_pfam
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Gianlucca G. Nicastro

#### 6. `PhyloBuilder` — Baixa

- **Assinatura:** `PhyloBuilder` — python, class
- **Local:** `lib/rotifer/devel/alpha/trsantos.py:539`–739 (201 linhas)
- **Última alteração (bruta):** 2025-10-13 · **excluindo ruído e importação em massa:** 2025-10-13
- **Commits que tocaram o corpo (1):** `f746b71`
- **Autores distintos (1):** Thiago Roberto dos Santos
- **Introdução (pickaxe):** [sem evidência]
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (16 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Thiago Roberto dos Santos`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Thiago Roberto dos Santos

#### 7. `get_seq_ids` — Baixa

- **Assinatura:** `get_seq_ids` — python, function
- **Local:** `share/rotifer/snakemake/rps/Snakefile.bkp:18`–18 (0 linhas)
- **Última alteração (bruta):** 2025-05-23 · **excluindo ruído e importação em massa:** 2025-05-23
- **Commits que tocaram o corpo (1):** `8619346`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_seq_ids("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - **outros: 15** — `share/rotifer/snakemake/neighborhood/Snakefile:11`, `share/rotifer/snakemake/neighborhood/Snakefile:29`, `share/rotifer/snakemake/rps/Snakefile.bkp:35`, `share/rotifer/snakemake/rps/Snakefile:18`, `share/rotifer/snakemake/rps/Snakefile:35`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp2:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp_working:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile2:19`…
- **Sinais encontrados:**
  - ARQUIVO NÃO PARSEIA: invalid syntax (linha 6)
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_seq_ids` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 8. `get_seq_ids` — Baixa

- **Assinatura:** `get_seq_ids` — python, function
- **Local:** `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp2:19`–19 (0 linhas)
- **Última alteração (bruta):** 2025-06-14 · **excluindo ruído e importação em massa:** 2025-06-14
- **Commits que tocaram o corpo (1):** `705301e`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_seq_ids("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - **outros: 15** — `share/rotifer/snakemake/neighborhood/Snakefile:11`, `share/rotifer/snakemake/neighborhood/Snakefile:29`, `share/rotifer/snakemake/rps/Snakefile.bkp:35`, `share/rotifer/snakemake/rps/Snakefile:18`, `share/rotifer/snakemake/rps/Snakefile:35`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp2:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp_working:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile2:19`…
- **Sinais encontrados:**
  - ARQUIVO NÃO PARSEIA: invalid syntax (linha 5)
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_seq_ids` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 9. `get_seq_ids` — Baixa

- **Assinatura:** `get_seq_ids` — python, function
- **Local:** `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp_working:19`–19 (0 linhas)
- **Última alteração (bruta):** 2026-02-26 · **excluindo ruído e importação em massa:** 2026-02-26
- **Commits que tocaram o corpo (1):** `558f52c`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_seq_ids("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - **outros: 15** — `share/rotifer/snakemake/neighborhood/Snakefile:11`, `share/rotifer/snakemake/neighborhood/Snakefile:29`, `share/rotifer/snakemake/rps/Snakefile.bkp:35`, `share/rotifer/snakemake/rps/Snakefile:18`, `share/rotifer/snakemake/rps/Snakefile:35`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp2:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile.bkp_working:37`, `share/rotifer/snakemake/rps_from_acc/Snakefile2:19`…
- **Sinais encontrados:**
  - ARQUIVO NÃO PARSEIA: invalid syntax (linha 5)
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_seq_ids` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 10. `reverse_domains(pid)` — Baixa

- **Assinatura:** `reverse_domains(pid)` (em `operon_fig.split_domain`) — python, method
- **Local:** `lib/rotifer/devel/alpha/draw.py:147`–157 (11 linhas)
- **Última alteração (bruta):** 2024-03-06 · **excluindo ruído e importação em massa:** 2024-03-06
- **Commits que tocaram o corpo (1):** `e5165d6`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def reverse_domains("`):** `e5165d6` — 2024-03-06 — rodolfoar · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/draw.py:176`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** rodolfoar (dono(s) do sandbox chamador)

#### 11. `fix_direction(df)` — Baixa

- **Assinatura:** `fix_direction(df)` (em `operon_fig.split_domain.reverse_domains`) — python, method
- **Local:** `lib/rotifer/devel/alpha/draw.py:148`–152 (5 linhas)
- **Última alteração (bruta):** 2024-03-06 · **excluindo ruído e importação em massa:** 2024-03-06
- **Commits que tocaram o corpo (1):** `e5165d6`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def fix_direction("`):** `e5165d6` — 2024-03-06 — rodolfoar · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/draw.py:153`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** rodolfoar (dono(s) do sandbox chamador)

#### 12. `rpsblast2table(psiblast_output, simple_profiledb)` — Baixa

- **Assinatura:** `rpsblast2table(psiblast_output, simple_profiledb)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1194`–1223 (30 linhas)
- **Última alteração (bruta):** 2024-09-04 · **excluindo ruído e importação em massa:** 2024-09-04
- **Commits que tocaram o corpo (2):** `f7ef7f6`, `0b74adf`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def rpsblast2table("`):** `0b74adf` — 2024-09-04 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:1439`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 13. `TMprediction(seqobj, predictior, cpu)` — Baixa

- **Assinatura:** `TMprediction(seqobj, predictior, cpu)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1249`–1275 (27 linhas)
- **Última alteração (bruta):** 2024-09-04 · **excluindo ruído e importação em massa:** 2024-09-04
- **Commits que tocaram o corpo (1):** `f7ef7f6`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def TMprediction("`):** `f7ef7f6` — 2024-09-04 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:1441`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 14. `get_plen(pidlist)` — Baixa

- **Assinatura:** `get_plen(pidlist)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1329`–1344 (16 linhas)
- **Última alteração (bruta):** 2024-09-10 · **excluindo ruído e importação em massa:** 2024-09-10
- **Commits que tocaram o corpo (2):** `a93375d`, `90dddc5`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_plen("`):** `90dddc5` — 2024-09-09 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/gian_func.py:1352`, `lib/rotifer/devel/alpha/gian_func.py:1354`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 15. `extract_organism_from_description(text)` — Baixa

- **Assinatura:** `extract_organism_from_description(text)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:1649`–1656 (8 linhas)
- **Última alteração (bruta):** 2024-09-10 · **excluindo ruído e importação em massa:** 2024-09-10
- **Commits que tocaram o corpo (1):** `a93375d`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def extract_organism_from_description("`):** `a93375d` — 2024-09-10 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:1464`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 16. `uniprot2fasta(uniprot_ids)` — Baixa

- **Assinatura:** `uniprot2fasta(uniprot_ids)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:3829`–3858 (30 linhas)
- **Última alteração (bruta):** 2025-03-14 · **excluindo ruído e importação em massa:** 2025-03-14
- **Commits que tocaram o corpo (1):** `12e6763`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def uniprot2fasta("`):** `12e6763` — 2025-03-14 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - url externa: https://rest.uniprot.org/uniprotkb/stream
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o comportamento desta função ainda é necessário, ou já foi absorvido por outro caminho de código?
- **Quem responde:** autor original: Gianlucca G. Nicastro

#### 17. `domain_seg(seqobj, ted_path, duplicates)` — Baixa

- **Assinatura:** `domain_seg(seqobj, ted_path, duplicates)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:3859`–3916 (58 linhas)
- **Última alteração (bruta):** 2025-03-14 · **excluindo ruído e importação em massa:** 2025-03-14
- **Commits que tocaram o corpo (1):** `12e6763`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def domain_seg("`):** `12e6763` — 2025-03-14 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:3921`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /netmnt/vast01/cbb01/proteinworld/People/gian/projects/TED/TED_database.db
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 18. `add_block_to_graph(graph, block_df, block_index, label_width, color_map, highlight_query, collapse_opposite_strand, font_size, left_pad, right_pad, spacer_width, show_row_label, seq_col)` — Baixa

- **Assinatura:** `add_block_to_graph(graph, block_df, block_index, label_width, color_map, highlight_query, collapse_opposite_strand, font_size, left_pad, right_pad, spacer_width, show_row_label, seq_col)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:1503`–1680 (178 linhas)
- **Última alteração (bruta):** 2026-09-09 · **excluindo ruído e importação em massa:** 2026-09-09
- **Commits que tocaram o corpo (11):** `4412cbd`, `57ee0c1`, `f71d83c`, `4e13e57`, `47d9d71`, `f848ed1`, `cfe26ed`, `7a525ef`…
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def add_block_to_graph("`):** `4412cbd` — 2026-06-23 — Eduardo Pereira Soares
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 4** — `lib/rotifer/devel/alpha/igem.py:104`, `lib/rotifer/devel/alpha/igem.py:1774`, `lib/rotifer/devel/alpha/igem.py:1866`, `lib/rotifer/devel/alpha/igem.py:2790`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (20 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** EdwardNSeven (dono(s) do sandbox chamador)

#### 19. `build_genome_overview_svg(extents, color_map, highlight_color, marker_color, track_width, left_margin, top_margin, row_height, track_height, font_size, label_lanes, label_lane_gap, max_labels_per_track, segment_length)` — Baixa

- **Assinatura:** `build_genome_overview_svg(extents, color_map, highlight_color, marker_color, track_width, left_margin, top_margin, row_height, track_height, font_size, label_lanes, label_lane_gap, max_labels_per_track, segment_length)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:2117`–2302 (186 linhas)
- **Última alteração (bruta):** 2026-09-04 · **excluindo ruído e importação em massa:** 2026-09-04
- **Commits que tocaram o corpo (5):** `5d9269a`, `0eac1f4`, `ed96cf2`, `0746b97`, `06b7aea`
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def build_genome_overview_svg("`):** `0eac1f4` — 2026-06-24 — EdwardNSeven
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 6** — `lib/rotifer/devel/alpha/igem.py:116`, `lib/rotifer/devel/alpha/igem.py:2336`, `lib/rotifer/devel/alpha/igem.py:2495`, `lib/rotifer/devel/alpha/igem.py:2499`, `lib/rotifer/devel/alpha/igem.py:2517`, `lib/rotifer/devel/alpha/igem.py:2537`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - url externa: http://www.w3.org/2000/svg
  - bloco comentado grande (13 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** EdwardNSeven (dono(s) do sandbox chamador)

#### 20. `build_genome_overview_interactive_html(extents, color_map, marker_color, highlight_color, base_track_height, segment_length)` — Baixa

- **Assinatura:** `build_genome_overview_interactive_html(extents, color_map, marker_color, highlight_color, base_track_height, segment_length)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:2317`–2475 (159 linhas)
- **Última alteração (bruta):** 2026-09-08 · **excluindo ruído e importação em massa:** 2026-09-08
- **Commits que tocaram o corpo (5):** `ed96cf2`, `5d9269a`, `9a47d03`, `0746b97`, `4e13e57`
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def build_genome_overview_interactive_html("`):** `ed96cf2` — 2026-06-24 — EdwardNSeven
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/igem.py:117`, `lib/rotifer/devel/alpha/igem.py:5709`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (14 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** EdwardNSeven (dono(s) do sandbox chamador)

#### 21. `build_scaled_block_svg(block_df, color_map, nucleotide_col, start_col, end_col, normalize_orientation, highlight_query, track_width, left_margin, right_margin, gene_height, font_size, min_gene_width, show_row_label, seq_col)` — Baixa

- **Assinatura:** `build_scaled_block_svg(block_df, color_map, nucleotide_col, start_col, end_col, normalize_orientation, highlight_query, track_width, left_margin, right_margin, gene_height, font_size, min_gene_width, show_row_label, seq_col)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:2927`–3200 (274 linhas)
- **Última alteração (bruta):** 2026-09-10 · **excluindo ruído e importação em massa:** 2026-09-10
- **Commits que tocaram o corpo (9):** `7959ce8`, `47d9d71`, `dc24cee`, `2bd7e3f`, `cfe26ed`, `7a525ef`, `522f287`, `06b7aea`…
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def build_scaled_block_svg("`):** `7959ce8` — 2026-08-27 — Eduardo Pereira Soares
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 5** — `lib/rotifer/devel/alpha/igem.py:121`, `lib/rotifer/devel/alpha/igem.py:1386`, `lib/rotifer/devel/alpha/igem.py:3206`, `lib/rotifer/devel/alpha/igem.py:3220`, `lib/rotifer/devel/alpha/igem.py:3230`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - url externa: http://www.w3.org/2000/svg
  - bloco comentado grande (23 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** EdwardNSeven (dono(s) do sandbox chamador)

#### 22. `build_html_report(df, output_file, title, group_col, org_col, label_col, rename_map, custom_colors, max_colors, ignore_domains, color_categories, nucleotide_col, start_col, end_col, length_col, seq_col, normalize_orientation, external_tools, table_include_seq, operon_kwargs, max_table_rows, work_dir, default_view, overview_segment_length, software_name, header_logo, footer_logo)` — Baixa

- **Assinatura:** `build_html_report(df, output_file, title, group_col, org_col, label_col, rename_map, custom_colors, max_colors, ignore_domains, color_categories, nucleotide_col, start_col, end_col, length_col, seq_col, normalize_orientation, external_tools, table_include_seq, operon_kwargs, max_table_rows, work_dir, default_view, overview_segment_length, software_name, header_logo, footer_logo)` — python, function
- **Local:** `lib/rotifer/devel/alpha/igem.py:5548`–5804 (257 linhas)
- **Última alteração (bruta):** 2026-09-10 · **excluindo ruído e importação em massa:** 2026-09-10
- **Commits que tocaram o corpo (17):** `ed96cf2`, `4e13e57`, `142da84`, `f28acfa`, `0eac1f4`, `522f287`, `5d9269a`, `3dc292d`…
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def build_html_report("`):** `0eac1f4` — 2026-06-24 — EdwardNSeven
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 15** — `lib/rotifer/devel/alpha/epsoares.py:1389`, `lib/rotifer/devel/alpha/igem.py:135`, `lib/rotifer/devel/alpha/igem.py:1791`, `lib/rotifer/devel/alpha/igem.py:1811`, `lib/rotifer/devel/alpha/igem.py:2138`, `lib/rotifer/devel/alpha/igem.py:216`, `lib/rotifer/devel/alpha/igem.py:2323`, `lib/rotifer/devel/alpha/igem.py:2369`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (20 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`EdwardNSeven`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Eduardo Pereira Soares; EdwardNSeven (dono(s) do sandbox chamador)

#### 23. `filter_fimo(fimoraw, gentab, repdist, gensize)` — Baixa

- **Assinatura:** `filter_fimo(fimoraw, gentab, repdist, gensize)` — python, function
- **Local:** `lib/rotifer/devel/alpha/malu.py:324`–419 (96 linhas)
- **Última alteração (bruta):** 2026-07-27 · **excluindo ruído e importação em massa:** 2026-07-27
- **Commits que tocaram o corpo (3):** `43e641a`, `cf1a588`, `ff25b90`
- **Autores distintos (1):** Maria Luiza Andreani
- **Introdução (pickaxe `git log -S "def filter_fimo("`):** `cf1a588` — 2026-07-20 — Maria Luiza Andreani
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/epsoares.py:1263`, `lib/rotifer/devel/alpha/epsoares.py:1331`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (14 linhas #)
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Maria Luiza Andreani`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Eduardo Pereira Soares (dono(s) do sandbox chamador)

#### 24. `calculate_positions(G, iterations)` — Baixa

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:111`–120 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `calculate_positions` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 25. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Baixa

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:122`–133 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `convert_to_cytoscape` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 26. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Baixa

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp:464`–496 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `get_network_community` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 27. `calculate_positions(G, iterations)` — Baixa

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:111`–120 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `calculate_positions` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 28. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Baixa

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:122`–133 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `convert_to_cytoscape` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 29. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Baixa

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.bkp2:641`–673 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `get_network_community` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 30. `calculate_positions(G, iterations)` — Baixa

- **Assinatura:** `calculate_positions(G, iterations)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:112`–121 (10 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def calculate_positions("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 8** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:760`, `lib/rotifer/devel/alpha/net_functions.bkp2:622`, `lib/rotifer/devel/alpha/net_functions.bkp:450`, `lib/rotifer/devel/alpha/net_functions.good.bkp:751`, `lib/rotifer/devel/alpha/net_functions.py:779`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:768`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:751`, `lib/rotifer/devel/alpha/net_functions2.py:751`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `calculate_positions` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 31. `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — Baixa

- **Assinatura:** `convert_to_cytoscape(G, positions, layout_key, scaling_factor, x_scale, y_scale)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:123`–134 (12 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def convert_to_cytoscape("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 14** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:630`, `lib/rotifer/devel/alpha/net_functions.bkp2:383`, `lib/rotifer/devel/alpha/net_functions.bkp2:545`, `lib/rotifer/devel/alpha/net_functions.bkp:376`, `lib/rotifer/devel/alpha/net_functions.good.bkp:441`, `lib/rotifer/devel/alpha/net_functions.good.bkp:642`, `lib/rotifer/devel/alpha/net_functions.py:455`, `lib/rotifer/devel/alpha/net_functions.py:656`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `convert_to_cytoscape` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 32. `get_network_community(G, community_to_color, weight, resolution_parameter)` — Baixa

- **Assinatura:** `get_network_community(G, community_to_color, weight, resolution_parameter)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.good.bkp:797`–829 (33 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_network_community("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 9** — `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:877`, `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py:972`, `lib/rotifer/devel/alpha/net_functions.bkp2:694`, `lib/rotifer/devel/alpha/net_functions.bkp:517`, `lib/rotifer/devel/alpha/net_functions.good.bkp:850`, `lib/rotifer/devel/alpha/net_functions.py:889`, `lib/rotifer/devel/alpha/net_functions.quebrada.py:884`, `lib/rotifer/devel/alpha/net_functions.working_20250627.py:850`…
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Nulo dentro do repositório: o arquivo é uma cópia de backup versionada, não importável por nome de módulo; o risco externo só existe se alguém importar o arquivo por caminho literal.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `get_network_community` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 33. `get_sheet(gc, name, index)` — Baixa

- **Assinatura:** `get_sheet(gc, name, index)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rfsouza.py:57`–62 (6 linhas)
- **Última alteração (bruta):** 2023-05-02 · **excluindo ruído e importação em massa:** 2023-05-02
- **Commits que tocaram o corpo (1):** `360a07f`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def get_sheet("`):** `360a07f` — 2023-05-01 — Robson Francisco de Souza
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/rfsouza.py:48`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Robson Francisco de Souza (dono(s) do sandbox chamador)

#### 34. `hmmscan_linear(sequences, file, models_path, cpus, columns, rename)` — Baixa

- **Assinatura:** `hmmscan_linear(sequences, file, models_path, cpus, columns, rename)` — python, function
- **Local:** `lib/rotifer/devel/alpha/epsoares.py:370`–431 (62 linhas)
- **Última alteração (bruta):** 2026-03-26 · **excluindo ruído e importação em massa:** 2026-03-26
- **Commits que tocaram o corpo (7):** `1d13cd4`, `829db2b`, `0221777`, `388df78`, `b9895b3`, `32f74dc`, `f67023f`
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def hmmscan_linear("`):** `388df78` — 2026-03-16 — EdwardNSeven · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /databases/pfam/Pfam-A.hmm
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/hmmer.py:160 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Eduardo Pereira Soares`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `hmmscan_linear` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 35. `annotate_seqobj(seqobj, ndf, cnt, selected_collumns)` — Baixa

- **Assinatura:** `annotate_seqobj(seqobj, ndf, cnt, selected_collumns)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:617`–652 (36 linhas)
- **Última alteração (bruta):** 2023-12-14 · **excluindo ruído e importação em massa:** 2023-12-14
- **Commits que tocaram o corpo (4):** `c2fb5ba`, `25663cc`, `db472c1`, `de526c6`
- **Autores distintos (2):** Gianlucca G. Nicastro, Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def annotate_seqobj("`):** `de526c6` — 2022-10-06 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `annotate_seqobj` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 36. `search2aln(df, coverage, evalue, id_with_coord, arch)` — Baixa

- **Assinatura:** `search2aln(df, coverage, evalue, id_with_coord, arch)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:726`–768 (43 linhas)
- **Última alteração (bruta):** 2023-09-28 · **excluindo ruído e importação em massa:** 2023-09-28
- **Commits que tocaram o corpo (3):** `de526c6`, `9066f37`, `6fad346`
- **Autores distintos (2):** Gianlucca G. Nicastro, Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def search2aln("`):** `de526c6` — 2022-10-06 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `search2aln` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 37. `full_annotate(seqobj, progress, batch_size, mirror, threads, after, before, eukaryotes, full)` — Baixa

- **Assinatura:** `full_annotate(seqobj, progress, batch_size, mirror, threads, after, before, eukaryotes, full)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:799`–820 (22 linhas)
- **Última alteração (bruta):** 2024-03-27 · **excluindo ruído e importação em massa:** 2024-03-27
- **Commits que tocaram o corpo (5):** `bb9f16e`, `25663cc`, `48d0c1e`, `a2e1c0d`, `27137a1`
- **Autores distintos (2):** Gianlucca G. Nicastro, Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def full_annotate("`):** `bb9f16e` — 2022-11-14 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `full_annotate` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 38. `get_alphafold_dssp(afid, sequence_object, id_to_guide, chain_id)` — Baixa

- **Assinatura:** `get_alphafold_dssp(afid, sequence_object, id_to_guide, chain_id)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:4063`–4116 (54 linhas)
- **Última alteração (bruta):** 2025-05-12 · **excluindo ruído e importação em massa:** 2025-05-12
- **Commits que tocaram o corpo (1):** `bfc33a9`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def get_alphafold_dssp("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/gian_func.py:4467`, `lib/rotifer/devel/alpha/gian_func.py:4469`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - url externa: https://alphafold.ebi.ac.uk/files/AF-{afid}-F1-model_v4.pdb
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado?
- **Quem responde:** Gianlucca G. Nicastro (dono(s) do sandbox chamador)

#### 39. `af_to_seq(seqobj)` — Baixa

- **Assinatura:** `af_to_seq(seqobj)` — python, function
- **Local:** `lib/rotifer/devel/alpha/idmapping_uniprot.py:215`–220 (6 linhas)
- **Última alteração (bruta):** 2023-08-02 · **excluindo ruído e importação em massa:** 2023-08-02
- **Commits que tocaram o corpo (1):** `922b31d`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def af_to_seq("`):** `922b31d` — 2023-08-02 — rodolfoar · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/db/uniprot/webapi/idmapping.py:530 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `af_to_seq` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 40. `network_dash(G, pos_dict, xx, function, color_schemes)` — Baixa

- **Assinatura:** `network_dash(G, pos_dict, xx, function, color_schemes)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions.py:137`–833 (697 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (3):** `1739875`, `ad06624`, `12e6763`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (26 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions2.py:137 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `network_dash` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 41. `network_dash(G, pos_dict, xx, function)` — Baixa

- **Assinatura:** `network_dash(G, pos_dict, xx, function)` — python, function
- **Local:** `lib/rotifer/devel/alpha/net_functions2.py:137`–794 (658 linhas)
- **Última alteração (bruta):** 2025-08-01 · **excluindo ruído e importação em massa:** 2025-08-01
- **Commits que tocaram o corpo (1):** `1739875`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def network_dash("`):** `ad06624` — 2025-01-23 — Gianlucca G. Nicastro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (25 linhas #)
- ⚠ **Nome ambíguo:** 8 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/net_functions.py:137 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 8 definições homônimas de `network_dash` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 42. `hmmer2aln(hmmout, df, threads, folder)` — Baixa

- **Assinatura:** `hmmer2aln(hmmout, df, threads, folder)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:147`–152 (6 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def hmmer2aln("`):** `463d209` — 2023-01-18 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:147 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `hmmer2aln` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 43. `hmmsearch_full2pandas(file, error_lines, keep_threshold)` — Baixa

- **Assinatura:** `hmmsearch_full2pandas(file, error_lines, keep_threshold)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:405`–452 (48 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def hmmsearch_full2pandas("`):** `db472c1` — 2022-10-17 — Gianlucca Gonçalves Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `hmmsearch_full2pandas` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 44. `hhr_to_aln(seqobj, hhr, database)` — Baixa

- **Assinatura:** `hhr_to_aln(seqobj, hhr, database)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:455`–548 (94 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def hhr_to_aln("`):** `df6d548` — 2022-08-02 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `hhr_to_aln` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 45. `annotate_seqobj(seqobj, df, cnt)` — Baixa

- **Assinatura:** `annotate_seqobj(seqobj, df, cnt)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:553`–580 (28 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def annotate_seqobj("`):** `de526c6` — 2022-10-06 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `annotate_seqobj` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 46. `alnxaln(seqobj, clustercol, minseq)` — Baixa

- **Assinatura:** `alnxaln(seqobj, clustercol, minseq)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:795`–832 (38 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def alnxaln("`):** `2dd9c9c` — 2022-11-29 — Gianlucca G. Nicastro · 3 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `alnxaln` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 47. `envelope_collection(df, alignment_method)` — Baixa

- **Assinatura:** `envelope_collection(df, alignment_method)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:865`–872 (8 linhas)
- **Última alteração (bruta):** 2023-02-08 · **excluindo ruído e importação em massa:** 2023-02-08
- **Commits que tocaram o corpo (1):** `5d72477`
- **Autores distintos (1):** Rodolfo Alvarenga Ribeiro
- **Introdução (pickaxe `git log -S "def envelope_collection("`):** `5d72477` — 2023-02-08 — Rodolfo Alvarenga Ribeiro · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:865 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `envelope_collection` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 48. `columns_to_positions(seqobj, start, end, model, evalue, complete, **kwargs)` — Baixa

- **Assinatura:** `columns_to_positions(seqobj, start, end, model, evalue, complete, **kwargs)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:1009`–1031 (23 linhas)
- **Última alteração (bruta):** 2023-07-13 · **excluindo ruído e importação em massa:** 2023-07-13
- **Commits que tocaram o corpo (3):** `0882fca`, `6a2a990`, `1967d75`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def columns_to_positions("`):** `0882fca` — 2023-07-11 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:1009 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `columns_to_positions` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 49. `count_arch(series, normalize, cut_off, count)` — Baixa

- **Assinatura:** `count_arch(series, normalize, cut_off, count)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:1041`–1063 (23 linhas)
- **Última alteração (bruta):** 2024-06-09 · **excluindo ruído e importação em massa:** 2024-06-09
- **Commits que tocaram o corpo (3):** `6a2a990`, `cfcd96f`, `6cd098c`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def count_arch("`):** `6a2a990` — 2023-07-13 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/beta/blast.py:1041 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `count_arch` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos

#### 50. `get_id_mapping_results_search(url)` — Baixa

- **Assinatura:** `get_id_mapping_results_search(url)` — python, function
- **Local:** `lib/rotifer/devel/alpha/uniprot.py:133`–157 (25 linhas)
- **Última alteração (bruta):** 2022-11-12 · **excluindo ruído e importação em massa:** 2022-11-12
- **Commits que tocaram o corpo (1):** `9a788b4`
- **Autores distintos (1):** Gianlucca Gonçalves Nicastro
- **Introdução (pickaxe `git log -S "def get_id_mapping_results_search("`):** `9a788b4` — 2022-11-12 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/idmapping_uniprot.py:246`, `lib/rotifer/devel/alpha/uniprot.py:179`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca Gonçalves Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_id_mapping_results_search` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca Gonçalves Nicastro; rodolfoar (dono(s) do sandbox chamador)_

#### 51. `af_to_seq(seqobj)` — Baixa

- **Assinatura:** `af_to_seq(seqobj)` — python, function
- **Local:** `lib/rotifer/db/uniprot/webapi/idmapping.py:530`–535 (6 linhas)
- **Última alteração (bruta):** 2023-08-02 · **excluindo ruído e importação em massa:** 2023-08-02
- **Commits que tocaram o corpo (1):** `602f66d`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def af_to_seq("`):** `922b31d` — 2023-08-02 — rodolfoar · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/idmapping_uniprot.py:215 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `af_to_seq` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** existe notebook, script pessoal ou pipeline fora do repositório que importe este símbolo? — _Prof. Robson de Souza / equipe do LEEP — só quem roda os workflows do laboratório pode responder_

#### 52. `html_highlight_consensus(s)` — Baixa

- **Assinatura:** `html_highlight_consensus(s)` — python, function
- **Local:** `lib/rotifer/devel/alpha/aln2pdf.py:379`–398 (20 linhas)
- **Última alteração (bruta):** 2025-11-18 · **excluindo ruído e importação em massa:** 2025-11-18
- **Commits que tocaram o corpo (2):** `bfc33a9`, `b8a044a`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def html_highlight_consensus("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 5** — `lib/rotifer/devel/alpha/aln2pdf.py:168`, `lib/rotifer/devel/alpha/aln2pdf.py:226`, `lib/rotifer/devel/alpha/aln2pdf.py:239`, `lib/rotifer/devel/alpha/gian_func.py:4424`, `lib/rotifer/devel/alpha/gian_func.py:4434`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: todo
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/gian_func.py:4349 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `html_highlight_consensus` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 53. `operon_fig(df, domain_dict)` — Baixa

- **Assinatura:** `operon_fig(df, domain_dict)` — python, function
- **Local:** `lib/rotifer/devel/alpha/draw.py:129`–223 (95 linhas)
- **Última alteração (bruta):** 2024-03-06 · **excluindo ruído e importação em massa:** 2024-03-06
- **Commits que tocaram o corpo (1):** `e5165d6`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def operon_fig("`):** `e5165d6` — 2024-03-06 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/igem.py:2`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/trsantos.py:777 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `operon_fig` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _EdwardNSeven (dono(s) do sandbox chamador)_

#### 54. `split_domain(df, column, fill, domain_dict, strand, after, before, remove_tm)` — Baixa

- **Assinatura:** `split_domain(df, column, fill, domain_dict, strand, after, before, remove_tm)` (em `operon_fig`) — python, method
- **Local:** `lib/rotifer/devel/alpha/draw.py:141`–193 (53 linhas)
- **Última alteração (bruta):** 2024-03-06 · **excluindo ruído e importação em massa:** 2024-03-06
- **Commits que tocaram o corpo (1):** `e5165d6`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def split_domain("`):** `e5165d6` — 2024-03-06 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 3** — `lib/rotifer/devel/alpha/draw.py:196`, `lib/rotifer/devel/alpha/gian_func.py:1917`, `lib/rotifer/devel/alpha/gian_func.py:2239`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/gian_func.py:1811 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `split_domain` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro; rodolfoar (dono(s) do sandbox chamador)_

#### 55. `padding_df(df, how)` — Baixa

- **Assinatura:** `padding_df(df, how)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:822`–836 (15 linhas)
- **Última alteração (bruta):** 2023-07-17 · **excluindo ruído e importação em massa:** 2023-07-17
- **Commits que tocaram o corpo (2):** `9d05e5c`, `8da21f0`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def padding_df("`):** `9d05e5c` — 2022-11-03 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:3504`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `padding_df` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 56. `html_highlight_consensus(s)` — Baixa

- **Assinatura:** `html_highlight_consensus(s)` — python, function
- **Local:** `lib/rotifer/devel/alpha/gian_func.py:4349`–4368 (20 linhas)
- **Última alteração (bruta):** 2025-05-12 · **excluindo ruído e importação em massa:** 2025-05-12
- **Commits que tocaram o corpo (1):** `bfc33a9`
- **Autores distintos (1):** Gianlucca G. Nicastro
- **Introdução (pickaxe `git log -S "def html_highlight_consensus("`):** `bfc33a9` — 2025-05-12 — Gianlucca G. Nicastro
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 5** — `lib/rotifer/devel/alpha/aln2pdf.py:168`, `lib/rotifer/devel/alpha/aln2pdf.py:226`, `lib/rotifer/devel/alpha/aln2pdf.py:239`, `lib/rotifer/devel/alpha/gian_func.py:4424`, `lib/rotifer/devel/alpha/gian_func.py:4434`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - marcadores: todo
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/aln2pdf.py:379 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Gianlucca G. Nicastro`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `html_highlight_consensus` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 57. `get_id_mapping_results_search(url)` — Baixa

- **Assinatura:** `get_id_mapping_results_search(url)` — python, function
- **Local:** `lib/rotifer/devel/alpha/idmapping_uniprot.py:134`–158 (25 linhas)
- **Última alteração (bruta):** 2023-02-09 · **excluindo ruído e importação em massa:** 2023-02-09
- **Commits que tocaram o corpo (1):** `7edc074`
- **Autores distintos (1):** Rodolfo Alvarenga Ribeiro
- **Introdução (pickaxe `git log -S "def get_id_mapping_results_search("`):** `9a788b4` — 2022-11-12 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/idmapping_uniprot.py:246`, `lib/rotifer/devel/alpha/uniprot.py:179`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`rodolfoar`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_id_mapping_results_search` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca Gonçalves Nicastro; rodolfoar (dono(s) do sandbox chamador)_

#### 58. `padding_df(df)` — Baixa

- **Assinatura:** `padding_df(df)` — python, function
- **Local:** `lib/rotifer/devel/alpha/rodolfo.py:782`–792 (11 linhas)
- **Última alteração (bruta):** 2023-01-18 · **excluindo ruído e importação em massa:** 2023-01-18
- **Commits que tocaram o corpo (1):** `463d209`
- **Autores distintos (1):** rodolfoar
- **Introdução (pickaxe `git log -S "def padding_df("`):** `9d05e5c` — 2022-11-03 — Gianlucca G. Nicastro · 4 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/gian_func.py:3504`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Robson Francisco de Souza`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `padding_df` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca G. Nicastro (dono(s) do sandbox chamador)_

#### 59. `operon_fig(df, group_col, label_col, org_col, output_file, max_colors, highlight_query, query_same_direction, font_size, ignore_domains)` — Baixa

- **Assinatura:** `operon_fig(df, group_col, label_col, org_col, output_file, max_colors, highlight_query, query_same_direction, font_size, ignore_domains)` — python, function
- **Local:** `lib/rotifer/devel/alpha/trsantos.py:777`–926 (150 linhas)
- **Última alteração (bruta):** 2026-06-03 · **excluindo ruído e importação em massa:** 2026-06-03
- **Commits que tocaram o corpo (1):** `e97f9d6`
- **Autores distintos (1):** Thiago Roberto dos Santos
- **Introdução (pickaxe `git log -S "def operon_fig("`):** `e5165d6` — 2024-03-06 — rodolfoar · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 1** — `lib/rotifer/devel/alpha/igem.py:2`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - bloco comentado grande (17 linhas #)
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/draw.py:129 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Baixo dentro do repositório, mas o arquivo é sandbox pessoal (`Thiago Roberto dos Santos`) e scripts do próprio dono fora do repo podem importá-lo — confirmar com ele antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `operon_fig` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _EdwardNSeven (dono(s) do sandbox chamador)_

#### 60. `hmmscan_linear(sequences, file, models_path, cpus, columns, rename)` — Baixa

- **Assinatura:** `hmmscan_linear(sequences, file, models_path, cpus, columns, rename)` — python, function
- **Local:** `lib/rotifer/devel/beta/hmmer.py:160`–221 (62 linhas)
- **Última alteração (bruta):** 2026-03-26 · **excluindo ruído e importação em massa:** 2026-03-26
- **Commits que tocaram o corpo (5):** `a1d6863`, `a791c58`, `9b8aeec`, `c9590ed`, `f67023f`
- **Autores distintos (2):** Eduardo Pereira Soares, EdwardNSeven
- **Introdução (pickaxe `git log -S "def hmmscan_linear("`):** `388df78` — 2026-03-16 — EdwardNSeven · 2 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - (c) devel/alpha: 0
  - (d) testes/doc: 0
  - outros: 0
- **Sinais encontrados:**
  - caminho absoluto: /databases/pfam/Pfam-A.hmm
- ⚠ **Nome ambíguo:** 2 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** `lib/rotifer/devel/alpha/epsoares.py:370 (homônimo — verificar se é a mesma função)`
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 2 definições homônimas de `hmmscan_linear` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** existe notebook, script pessoal ou pipeline fora do repositório que importe este símbolo? — _Prof. Robson de Souza / equipe do LEEP — só quem roda os workflows do laboratório pode responder_

#### 61. `get_id_mapping_results_search(url)` — Baixa

- **Assinatura:** `get_id_mapping_results_search(url)` — python, function
- **Local:** `lib/rotifer/db/uniprot/webapi/idmapping.py:461`–485 (25 linhas)
- **Última alteração (bruta):** 2023-08-02 · **excluindo ruído e importação em massa:** 2023-08-02
- **Commits que tocaram o corpo (1):** `602f66d`
- **Autores distintos (1):** Robson Francisco de Souza
- **Introdução (pickaxe `git log -S "def get_id_mapping_results_search("`):** `9a788b4` — 2022-11-12 — Gianlucca Gonçalves Nicastro · 6 commits mexeram nessa assinatura
- **Call sites:**
  - (a) núcleo: 0
  - (b) CLI: 0
  - **(c) devel/alpha: 2** — `lib/rotifer/devel/alpha/idmapping_uniprot.py:246`, `lib/rotifer/devel/alpha/uniprot.py:179`
  - (d) testes/doc: 0
  - outros: 0
- **Sinais textuais de obsolescência:** nenhum
- ⚠ **Nome ambíguo:** 3 definições homônimas no repositório em arquivos distintos; as contagens acima podem pertencer a outra definição
- **Substituta sugerida:** [sem evidência]
- **Risco de quebra retroativa:** Módulo do núcleo importável: mesmo sem call site interno, um notebook ou script pessoal fora do repo pode importar este símbolo — exige confirmação humana antes de remover.
- **Pergunta a responder antes de qualquer ação:** as 3 definições homônimas de `get_id_mapping_results_search` são a mesma função copiada ou funções diferentes? Se forem diferentes, a qual delas pertencem os call sites contados acima?
- **Quem responde:** quem escreveu o código — as contagens desta ficha não distinguem homônimos
- **Pergunta adicional:** o sandbox que a consome ainda está em uso, ou é scratchpad abandonado? — _Gianlucca Gonçalves Nicastro; rodolfoar (dono(s) do sandbox chamador)_
- **Pergunta adicional:** existe notebook, script pessoal ou pipeline fora do repositório que importe este símbolo? — _Prof. Robson de Souza / equipe do LEEP — só quem roda os workflows do laboratório pode responder_

---

## Núcleo consumido só por devel/alpha

Funções definidas **fora** de `lib/rotifer/devel/alpha/` cujos **únicos** consumidores internos estão dentro do sandbox pessoal. Não são candidatas a remoção pelo critério do eixo 2 — são um alerta de acoplamento: código do núcleo mantido vivo apenas por scratchpads individuais.

| Função                          | Definida em                                      | Call sites em devel/alpha                                                                                                                                                                                                                         | Dono aparente do sandbox chamador (`git log --format=%an`)                                                                                 |
| ------------------------------- | ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `update_database`               | `lib/rotifer/db/local/ete3.py:33`                | `lib/rotifer/devel/alpha/trsantos.py:302`, `lib/rotifer/devel/alpha/trsantos.py:475`                                                                                                                                                              | `trsantos.py` → Thiago Roberto dos Santos (19 commits)                                                                                     |
| `get_id_mapping_results_search` | `lib/rotifer/db/uniprot/webapi/idmapping.py:461` | `lib/rotifer/devel/alpha/idmapping_uniprot.py:246`, `lib/rotifer/devel/alpha/uniprot.py:179`                                                                                                                                                      | `idmapping_uniprot.py` → rodolfoar (9 commits); `uniprot.py` → Gianlucca Gonçalves Nicastro (1 commits)                                    |
| `padding_df`                    | `lib/rotifer/devel/beta/blast.py:782`            | `lib/rotifer/devel/alpha/gian_func.py:3504`                                                                                                                                                                                                       | `gian_func.py` → Gianlucca G. Nicastro (119 commits)                                                                                       |
| `alnclu`                        | `lib/rotifer/devel/beta/blast.py:952`            | `lib/rotifer/devel/alpha/rodolfo.py:959`                                                                                                                                                                                                          | `rodolfo.py` → Robson Francisco de Souza (51 commits)                                                                                      |
| `trim`                          | `lib/rotifer/devel/beta/sequence.py:1595`        | `lib/rotifer/devel/alpha/gian_func.py:948`, `lib/rotifer/devel/alpha/igem.py:4878`, `lib/rotifer/devel/alpha/igem.py:4892`, `lib/rotifer/devel/alpha/igem.py:5067`, `lib/rotifer/devel/alpha/igem.py:5107`, `lib/rotifer/devel/alpha/malu.py:637` | `gian_func.py` → Gianlucca G. Nicastro (119 commits); `igem.py` → EdwardNSeven (43 commits); `malu.py` → Maria Luiza Andreani (22 commits) |
| `select_neighbors`              | `lib/rotifer/genome/data.py:355`                 | `lib/rotifer/devel/alpha/draw.py:159`, `lib/rotifer/devel/alpha/epsoares.py:153`                                                                                                                                                                  | `draw.py` → rodolfoar (8 commits); `epsoares.py` → Eduardo Pereira Soares (113 commits)                                                    |

Observações por caso:

- **`update_database`** (`lib/rotifer/db/local/ete3.py:33`) — última alteração 2022-12-20.
- **`get_id_mapping_results_search`** (`lib/rotifer/db/uniprot/webapi/idmapping.py:461`) — última alteração 2023-08-02. ⚠ nome ambíguo (3 definições homônimas) — confirmar que o call site resolve para esta definição antes de agir.
- **`padding_df`** (`lib/rotifer/devel/beta/blast.py:782`) — última alteração 2025-07-14. ⚠ nome ambíguo (3 definições homônimas) — confirmar que o call site resolve para esta definição antes de agir.
- **`alnclu`** (`lib/rotifer/devel/beta/blast.py:952`) — última alteração 2025-07-14. ⚠ nome ambíguo (2 definições homônimas) — confirmar que o call site resolve para esta definição antes de agir.
- **`trim`** (`lib/rotifer/devel/beta/sequence.py:1595`) — última alteração 2025-08-01.
- **`select_neighbors`** (`lib/rotifer/genome/data.py:355`) — última alteração 2019-12-11.

**Donos dos sandboxes em `lib/rotifer/devel/alpha/`** (autor com mais commits em cada arquivo):

| Arquivo                                                                 | Dono aparente                | Commits | Último commit |
| ----------------------------------------------------------------------- | ---------------------------- | ------- | ------------- |
| `lib/rotifer/devel/alpha/__init__.py`                                   | Gianlucca Gonçalves Nicastro | 9       | 2022-06-06    |
| `lib/rotifer/devel/alpha/aln2pdf.py`                                    | Gianlucca G. Nicastro        | 3       | 2025-11-18    |
| `lib/rotifer/devel/alpha/collection.py`                                 | Robson Francisco de Souza    | 5       | 2023-10-11    |
| `lib/rotifer/devel/alpha/draw.py`                                       | rodolfoar                    | 8       | 2025-04-03    |
| `lib/rotifer/devel/alpha/epsoares.py`                                   | Eduardo Pereira Soares       | 113     | 2026-09-09    |
| `lib/rotifer/devel/alpha/genome/draw.py`                                | Robson Francisco de Souza    | 1       | 2024-05-02    |
| `lib/rotifer/devel/alpha/gian_func.py`                                  | Gianlucca G. Nicastro        | 119     | 2026-04-13    |
| `lib/rotifer/devel/alpha/idmapping_uniprot.py`                          | rodolfoar                    | 9       | 2024-05-29    |
| `lib/rotifer/devel/alpha/igem.py`                                       | EdwardNSeven                 | 43      | 2026-09-09    |
| `lib/rotifer/devel/alpha/malu.py`                                       | Maria Luiza Andreani         | 22      | 2026-09-06    |
| `lib/rotifer/devel/alpha/mapper.py`                                     | Robson Francisco de Souza    | 1       | 2022-06-22    |
| `lib/rotifer/devel/alpha/net_functions.backup_meioquebrada_20250708.py` | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.bkp`                             | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.bkp2`                            | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.good.bkp`                        | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.py`                              | Gianlucca G. Nicastro        | 3       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.quebrada.py`                     | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions.working_20250627.py`             | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/net_functions2.py`                             | Gianlucca G. Nicastro        | 1       | 2025-08-01    |
| `lib/rotifer/devel/alpha/rfsouza.py`                                    | Robson Francisco de Souza    | 1       | 2023-05-01    |
| `lib/rotifer/devel/alpha/rodolfo.py`                                    | Robson Francisco de Souza    | 51      | 2025-03-11    |
| `lib/rotifer/devel/alpha/sequence/__init__.py`                          | rodolfoar                    | 1       | 2024-02-29    |
| `lib/rotifer/devel/alpha/sequence/complete.py`                          | rodolfoar                    | 2       | 2024-03-01    |
| `lib/rotifer/devel/alpha/snakemake.py`                                  | Gianlucca G. Nicastro        | 6       | 2026-04-13    |
| `lib/rotifer/devel/alpha/snakemake2.py`                                 | Gianlucca G. Nicastro        | 1       | 2025-05-12    |
| `lib/rotifer/devel/alpha/trsantos.py`                                   | Thiago Roberto dos Santos    | 19      | 2026-06-03    |
| `lib/rotifer/devel/alpha/uniprot.py`                                    | Gianlucca Gonçalves Nicastro | 1       | 2022-11-12    |

---

## Falsos positivos descartados

As 37 definições abaixo passaram por todos os filtros mecânicos — zero call site em núcleo, CLI, testes e documentação — mas **não** são código legado. Cada uma traz a evidência que a resgata.

### função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada — 14

| Função                  | Local                                       | Por que não é legado                                                                                                                       |
| ----------------------- | ------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `acc2pfam`              | `etc/profile.d/acc2pfam.zsh:1`              | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `acc2profiledb`         | `etc/profile.d/acc2profiledb.sh:1`          | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `acc2profiledb`         | `etc/profile.d/acc2profiledb.zsh:1`         | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `aln_order_by_tree`     | `etc/profile.d/aln_order_by_tree.sh:1`      | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `aln_order_by_tree`     | `etc/profile.d/aln_order_by_tree.zsh:1`     | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `distinct_dom_per_prot` | `etc/profile.d/distinct_dom_per_prot.sh:1`  | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `distinct_dom_per_prot` | `etc/profile.d/distinct_dom_per_prot.zsh:1` | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `dom_count`             | `etc/profile.d/dom_count.sh:1`              | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `dom_count`             | `etc/profile.d/dom_count.zsh:1`             | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `dom_per_prot`          | `etc/profile.d/dom_per_prot.sh:1`           | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `dom_per_prot`          | `etc/profile.d/dom_per_prot.zsh:1`          | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `runhmmscan`            | `etc/profile.d/runhmmscan.sh:1`             | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `tdesc`                 | `etc/profile.d/tdesc.sh:1`                  | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |
| `tdesc`                 | `etc/profile.d/tdesc.zsh:1`                 | função de shell publicada em etc/profile.d/ — é API interativa carregada no shell do pesquisador; ausência de call site interno é esperada |

### arquivo usa dispatch dinâmico — 12

| Função                        | Local                                      | Por que não é legado                                                                                         |
| ----------------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| `ipgs_to_dict`                | `lib/rotifer/db/methods.py:183`            | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `drop_release`                | `lib/rotifer/db/uniprot/clickhouse.py:443` | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `add_hhpred`                  | `lib/rotifer/devel/beta/sequence.py:730`   | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `to_bioalign`                 | `lib/rotifer/devel/beta/sequence.py:1060`  | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `add_jpred`                   | `lib/rotifer/devel/beta/sequence.py:1624`  | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `_to_df_styleTEX`             | `lib/rotifer/devel/beta/sequence.py:2063`  | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `edit`                        | `lib/rotifer/devel/beta/sequence.py:2208`  | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `fetch_neighbors`             | `lib/rotifer/devel/beta/sequence.py:2261`  | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `series_to_compact_frequency` | `lib/rotifer/genome/data.py:338`           | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `select_neighbors`            | `lib/rotifer/genome/data.py:355`           | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `full_neighborhood`           | `lib/rotifer/genome/data.py:541`           | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |
| `jaccard`                     | `lib/rotifer/genome/data.py:890`           | arquivo usa dispatch dinâmico (globals()/getattr/eval com string) — call site não é rastreável estaticamente |

### hook de framework — 3

| Função         | Local                              | Por que não é legado                                                                                   |
| -------------- | ---------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `emit`         | `lib/rotifer/core/logger.py:12`    | hook de framework (emit é invocado internamente por pandas/logging, nunca por nome no código)          |
| `_constructor` | `lib/rotifer/genome/data.py:73`    | hook de framework (\_constructor é invocado internamente por pandas/logging, nunca por nome no código) |
| `_constructor` | `lib/rotifer/seq/alignment.py:196` | hook de framework (\_constructor é invocado internamente por pandas/logging, nunca por nome no código) |

### callback de framework Perl/DBIx::Class — 3

| Função         | Local                                                                       | Por que não é legado                                                                                                  |
| -------------- | --------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `DESTROY`      | `perl/lib/Rotifer/DB/NCBI/Taxonomy.pm:996`                                  | callback de framework Perl/DBIx::Class (DESTROY é chamado implicitamente pelo interpretador ou pela classe base)      |
| `can_call_new` | `perl/lib/Rotifer/DBIC/AnnotationDB/Component/BiosequenceComponent.pm:208`  | callback de framework Perl/DBIx::Class (can_call_new é chamado implicitamente pelo interpretador ou pela classe base) |
| `store_column` | `perl/lib/Rotifer/DBIC/AnnotationDB/Component/DynamicSubclassByTerm.pm:161` | callback de framework Perl/DBIx::Class (store_column é chamado implicitamente pelo interpretador ou pela classe base) |

### invocada por dispatch dinâmico no próprio arquivo — 3

| Função        | Local              | Por que não é legado                                                       |
| ------------- | ------------------ | -------------------------------------------------------------------------- |
| `runhmmscan`  | `bin/scanseqs:186` | invocada por dispatch dinâmico no próprio arquivo (construção "run${...}") |
| `runrpsblast` | `bin/scanseqs:224` | invocada por dispatch dinâmico no próprio arquivo (construção "run${...}") |
| `runpsiblast` | `bin/scanseqs:243` | invocada por dispatch dinâmico no próprio arquivo (construção "run${...}") |

### decorado como property/setter — acessado como atributo, não como chamada — 2

| Função            | Local                                     | Por que não é legado                                                     |
| ----------------- | ----------------------------------------- | ------------------------------------------------------------------------ |
| `freq_table`      | `lib/rotifer/devel/beta/sequence.py:989`  | decorado como property/setter — acessado como atributo, não como chamada |
| `from_seqrecords` | `lib/rotifer/devel/beta/sequence.py:2322` | decorado como property/setter — acessado como atributo, não como chamada |

### Classes inteiras de falso positivo filtradas antes da pontuação

| Classe                         | Quantidade | Evidência                                                                                                             |
| ------------------------------ | ---------- | --------------------------------------------------------------------------------------------------------------------- |
| Métodos `__dunder__`           | 171        | Invocados pelo interpretador via protocolo (`__init__`, `__getitem__`, `__enter__`…), nunca por nome no código-fonte. |
| Funções `pytest` em `test/`    | 413        | Coletadas por convenção de nome pelo runner; ausência de call site é o funcionamento normal do framework.             |
| Exemplos e templates em `doc/` | 1          | `doc/examples/basic_rotifer_script_template.py` e similares existem para ser copiados, não chamados.                  |

### Casos verificados individualmente que merecem registro

**`runhmmscan`, `runrpsblast`, `runpsiblast`, `runmmseqs`, `runhmmsearch`, `runphobius` — `bin/scanseqs`.** Nenhum tem call site textual. Evidência que os resgata: `bin/scanseqs:309` monta o nome da função em tempo de execução —

```bash
run="run${pipeline[1]}"
```

— onde `pipeline[1]` vem das especificações `nome=programa:banco` das linhas 11–13 (`aravind=rpsblast:…`, `pfam=hmmscan:…`). Todo o conjunto `run*` é despachado dinamicamente. Nenhuma varredura estática de call site poderia encontrá-los.

**`etc/profile.d/*.sh` e `etc/profile.d/*.zsh` (14 funções: `acc2pfam`, `acc2profiledb`, `aln_order_by_tree`, `distinct_dom_per_prot`, `dom_count`, `dom_per_prot`, `tdesc`, `runhmmscan`).** São funções de shell publicadas em `profile.d/` — o diretório existe justamente para ser carregado no shell interativo do pesquisador. A ausência de call site interno é o comportamento esperado, e é exatamente o caso que a restrição de compatibilidade retroativa do projeto protege. Cada uma existe em duas cópias (`.sh` e `.zsh`), o que é redundância deliberada de dialeto, não duplicação legada.

**`DESTROY` (`perl/lib/Rotifer/DB/NCBI/Taxonomy.pm:996`).** Método especial do Perl: chamado pelo interpretador na destruição do objeto. Remover causaria vazamento de recurso, não limpeza.

**`can_call_new` (`.../BiosequenceComponent.pm:208`) e `store_column` (`.../DynamicSubclassByTerm.pm:161`).** São _overrides_ da API de componentes do DBIx::Class, invocados pela classe base via `load_components`. O chamador está na dependência, não no repositório.

**`_constructor` (`lib/rotifer/genome/data.py:73`, `lib/rotifer/seq/alignment.py:196`) e `emit` (`lib/rotifer/core/logger.py:12`).** `_constructor` é a propriedade que o pandas consulta internamente ao propagar subclasses de `DataFrame`; `emit` é o método que `logging.Handler` chama. Ambos são contratos de framework.

---

## Anexo — o que esta auditoria _não_ prova

1. **"Sem referência interna" não é "sem uso".** O repositório não contém os scripts pessoais, notebooks e pipelines que os pesquisadores do LEEP rodam fora dele. Toda função do núcleo (`lib/`) e todo módulo Perl instalável listado aqui pode ter consumidores invisíveis a esta varredura. Por isso a camada **Alta** é dominada por cópias de backup, arquivos que não compilam e código de sandbox — casos em que o consumo externo é fisicamente improvável — e não por "função antiga do núcleo sem chamador".
2. **A resolução de nomes é textual.** Não há análise de escopo, herança ou import. Homônimos inflam a contagem de uso (direção segura) e ambiguidade foi sempre sinalizada.
3. **`git blame` atribui a linha ao último commit que a tocou.** A "data de introdução" derivada do blame é a da linha sobrevivente mais antiga; a introdução real vem do pickaxe, reportado separadamente em cada ficha. Quando as duas divergem, a divergência costuma significar que a função foi reescrita.
4. **Nenhuma execução.** Não foram instaladas dependências, não foi importado nenhum módulo e nenhum teste foi rodado. Os cinco arquivos que não compilam foram detectados por `ast.parse`, não por execução.
