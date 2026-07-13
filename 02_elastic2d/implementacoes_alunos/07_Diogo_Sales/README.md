# Pressão Distribuída no Contorno e Cálculo de Tensões (Elementos T3/Q4)

**Grupo 07 — Diogo Sales e Pedro López**

---

## 1. Descrição

Esta implementação estende o programa base de MEF elástico 2D (elementos T3 e Q4, estado plano de tensão/deformação) com duas funcionalidades que não existiam no código de referência: **aplicação de pressão distribuída no contorno** e **cálculo/exportação de tensões por elemento**.

O programa original só aceitava forças concentradas aplicadas diretamente nos nós. Cargas reais, no entanto, costumam ser distribuídas ao longo de uma superfície ou aresta, exemplo, pressão de um fluido. Para incorporar esse tipo de carregamento, foi necessário converter a pressão, informada ao longo de uma aresta do contorno, em um par de forças nodais **equivalentes**, aplicadas nos dois nós daquela aresta.

Além disso, o programa de referência resolvia o sistema $[K]\{U\}=\{F\}$ e escrevia os deslocamentos, mas não calculava as tensões resultantes nem exportava nenhum resultado por elemento para visualização gráfica. Foi implementada uma rotina de pós-processamento que recupera a tensão ($\sigma_x, \sigma_y, \tau_{xy}, \sigma_z$) em cada elemento a partir do campo de deslocamentos resolvido, e o resultado foi exportado no formato VTK (`CELL_DATA`) para visualização direta no ParaView.

Nenhuma rotina do núcleo do programa de referência (elementos `elmt01`–`elmt04`, montagem `pform`/`addstf`, solver PCG) foi alterada — a implementação estende o programa por fora, reaproveitando essas peças já validadas.

---

## 2. Formulação teórica

### 2.1 Forças nodais equivalentes a uma pressão de borda

Para uma aresta linear de 2 nós $(n_i, n_j)$, com pressão uniforme $p$, a força nodal equivalente em cada nó é dada pelo princípio dos trabalhos virtuais:

$$F_a = \int_{\Gamma} N_a(\xi)\cdot t \; d\Gamma, \qquad a = 1,2$$

onde $N_1(\xi) = \dfrac{1-\xi}{2}$ e $N_2(\xi) = \dfrac{1+\xi}{2}$ são as funções de forma lineares da aresta ($\xi \in [-1,1]$), e $t = p\cdot n$ é o vetor de tração — a pressão escalar multiplicada pela normal unitária que aponta para fora do domínio:

$$n = \frac{(dy,\, -dx)}{L}, \qquad (dx,dy) = (x_j,y_j) - (x_i,y_i)$$

A normal só aponta corretamente para fora se os nós forem informados na ordem **anti-horária** do contorno (convenção: o sólido fica sempre à esquerda de quem caminha de $n_i$ para $n_j$).

A integral é resolvida numericamente por **quadratura de Gauss de 2 pontos**:

$$\xi_{1,2} = \mp\frac{1}{\sqrt3}, \qquad w_{1,2} = 1$$

Como o integrando ($N_a \cdot t$) é de grau 1, a quadratura de 2 pontos é **exata** — o resultado coincide com a fórmula fechada clássica $F = p\cdot L/2$ em cada nó, mas obtido de forma explícita a partir das funções de forma, o que permite generalizar facilmente para pressão variável ao longo da aresta (bastaria interpolar $p(\xi) = N_1 p_1 + N_2 p_2$).

### 2.2 Cálculo de tensões

A partir do campo de deslocamentos $u$ resolvido pelo solver, a deformação em cada elemento é obtida por:

$$\varepsilon_{xx} = \frac{\partial u_x}{\partial x}, \qquad
\varepsilon_{yy} = \frac{\partial u_y}{\partial y}, \qquad
\gamma_{xy} = \frac{\partial u_y}{\partial x} + \frac{\partial u_x}{\partial y}$$

calculadas via as derivadas das funções de forma do próprio elemento (matriz $B$: $\varepsilon = B\, u_e$). Para o T3 (deformação constante), essas derivadas são avaliadas uma única vez por elemento; para o Q4, são avaliadas no centro do elemento ($r=s=0$, equivalente a um único ponto de Gauss), produzindo um valor representativo de tensão por elemento.

A tensão é obtida pela lei de Hooke generalizada, $\sigma = D\varepsilon$, com a matriz constitutiva $D$ dependendo do estado plano considerado:

---

## 3. Como compilar e rodar

### 3.1 Pré-requisitos
- `gfortran` (compilador Fortran, suporta a sintaxe legada usada no arquivo)

### 3.2 Compilação
```bash
gfortran -o mef Mef-Entrega.f
```

### 3.3 Execução

O programa pede três nomes de arquivo, em ordem: dados de entrada, saída em texto, saída em VTK.

**Interativo:**
```bash
./mef
```
```
Arquivo de dados: caso1_uniaxial.dat
Arquivo de saida: saida1.txt
Arquivo VTK de saida (ex: resultado.vtk): caso1.vtk
```

### 3.3 Formato do bloco de pressão no arquivo `.dat`

Após as seções de condições de contorno e forças concentradas (já existentes no formato original), o novo bloco de pressão segue o formato:
```
ni   nj   p
```
uma linha por aresta carregada, terminando com o sentinela `0 0 0.0`. Os nós devem ser informados na ordem anti-horária do contorno.

---

## 4. Exemplo

Caso de teste `caso1_uniaxial.dat` — patch test uniaxial:

- Malha retangular 2,0 × 1,0 (T3, Estado Plano de Tensão), 15 nós, 16 elementos.
- Apoio tipo rolete na face esquerda (nós 1, 6, 11 com $u_x=0$; nó 1 também com $u_y=0$ para eliminar movimento de corpo rígido, sem restringir o efeito de Poisson).
- Pressão $p=100$ aplicada na face direita (arestas 5→10 e 10→15).

**Trecho relevante do `.dat`:**
```
 5   10   100.00
10   15   100.00
 0    0     0.0
```

**Resultado obtido** (idêntico em todos os 16 elementos):

| Campo | Esperado | Obtido |
|---|---|---|
| $\sigma_x$ | 100,0 | 100,0 (exato) |
| $\sigma_y$ | 0 | ≈ 0 (ruído numérico ~1e-6) |
| $\tau_{xy}$ | 0 | ≈ 0 (ruído numérico ~1e-6) |

O resultado é visualizável abrindo o `caso1.vtk` gerado no ParaView e colorindo por `sigma_x` (`CELL_DATA`) — a peça inteira aparece com uma única cor, confirmando o campo de tensão uniforme.

---

## 5. Verificação

A implementação foi validada em três níveis crescentes de exigência, além de dois testes de robustez de entrada:

1. **Patch test uniaxial** (`caso1_uniaxial.dat`) — pressão em uma única face. Resultado: $\sigma_x=100{,}0$ exato em todos os elementos, $\sigma_y\approx\tau_{xy}\approx 0$. Valida a conversão básica de pressão em força nodal equivalente.

2. **Patch test biaxial** (`caso2_biaxial.dat`) — pressão simultânea em duas faces perpendiculares. Resultado: $\sigma_x=\sigma_y=100{,}0$ exato, $\tau_{xy}\approx 0$. Valida a superposição de contribuições de arestas diferentes no mesmo nó (acumulação `+=` em `edgfor`).

Após cada alteração incremental do código, os resultados foram comparados com a versão anterior para confirmar que permaneciam idênticos (mesma norma de energia, mesmas tensões) antes de prosseguir — prática de não regressão.

---

## 6. Referências

- Notas de Aula do Prof. Fernando L. B. Ribeiro - INTRODUÇÃO AO MÉTODO DOS ELEMENTOS FINITOS
- Documentação do formato VTK legado (ASCII), para a estrutura dos blocos `POINTS`, `CELLS`, `CELL_TYPES`, `POINT_DATA` e `CELL_DATA`.
