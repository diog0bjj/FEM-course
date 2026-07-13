# Pressão Distribuída no Contorno e Cálculo de Tensões (Elementos T3/Q4)

**Grupo 07 — Diogo Sales e Pedro López**

---

## 1. Descrição

Este trabalho estende o programa base de MEF elástico 2D (elementos T3 e Q4, estado plano de tensão e de deformação) com dois recursos que não existiam no código de referência: aplicação de **pressão distribuída no contorno** e **cálculo e exportação das tensões por elemento**.

O programa original só aceitava forças concentradas nos nós. Como cargas reais costumam ser distribuídas ao longo de uma aresta (a pressão de um fluido, por exemplo), foi preciso converter essa pressão em forças nodais **equivalentes**, aplicadas nos dois nós da aresta carregada.

Além disso, o código de referência resolvia $[K]\{U\}=\{F\}$ e escrevia apenas os deslocamentos — não calculava as tensões nem exportava resultados por elemento. Implementamos uma rotina de pós-processamento que recupera a tensão em cada elemento a partir do campo de deslocamentos, com saída em coordenadas cartesianas ($\sigma_x, \sigma_y, \tau_{xy}, \sigma_z$) ou polares ($\sigma_r, \sigma_\theta, \tau_{r\theta}, \sigma_z$), escolhidas por uma flag no arquivo de entrada. O resultado é gravado em formato VTK (`CELL_DATA`) para visualização direta no ParaView.

As rotinas de cálculo do núcleo (elementos `elmt01`–`elmt04`, montagem `pform`/`addstf`, solver PCG) não foram alteradas. As mudanças ficaram nas rotinas de leitura e saída (`rdata`, `wvtk`), estendidas para os novos recursos, e nas rotinas novas de carga distribuída e de recuperação de tensões.

---

## 2. Formulação

### 2.1 Forças nodais equivalentes a uma pressão de borda

Para uma aresta linear de 2 nós $(n_i, n_j)$ com pressão uniforme $p$, a força equivalente em cada nó vem do princípio dos trabalhos virtuais:

$$F_a = \int_{\Gamma} N_a(\xi)\, t \; d\Gamma, \qquad a = 1,2$$

com $N_1(\xi) = \dfrac{1-\xi}{2}$, $N_2(\xi) = \dfrac{1+\xi}{2}$ e $t = p\, n$, onde $n$ é a normal unitária externa à aresta:

$$n = \frac{(dy,\, -dx)}{L}, \qquad (dx,dy) = (x_j,y_j) - (x_i,y_i)$$

A tração $t = p\,n$ é um vetor, então a força em cada nó também tem componentes em $x$ e $y$. A normal só aponta para fora se os nós forem informados na ordem **anti-horária** do contorno (o sólido fica sempre à esquerda de quem caminha de $n_i$ para $n_j$). Uma pressão que comprime a face entra com sinal negativo.

A integral é resolvida por **quadratura de Gauss de 2 pontos** ($\xi = \mp 1/\sqrt3$, pesos iguais a 1). Como o integrando ($N_a \cdot t$) é de grau 1, a quadratura é exata e reproduz o resultado clássico $F = p\,L/2$ em cada nó — mas obtido a partir das funções de forma, o que deixa a rotina pronta para o caso de pressão variável ao longo da aresta.

### 2.2 Cálculo de tensões

A deformação em cada elemento é obtida das derivadas das funções de forma aplicadas aos deslocamentos nodais ($\varepsilon = B\, u_e$):

$$\varepsilon_{xx} = \frac{\partial u_x}{\partial x}, \qquad
\varepsilon_{yy} = \frac{\partial u_y}{\partial y}, \qquad
\gamma_{xy} = \frac{\partial u_y}{\partial x} + \frac{\partial u_x}{\partial y}$$

No T3 (deformação constante) essas derivadas são avaliadas uma única vez por elemento. No Q4 elas variam dentro do elemento, então são avaliadas no centro ($r=s=0$, equivalente a um ponto de Gauss), produzindo um valor representativo de tensão por elemento.

A tensão vem da lei de Hooke, $\sigma = D\,\varepsilon$, com a matriz $D$ dependendo do estado plano:

**Estado plano de tensão (EPT):**

$$D = \frac{E}{1-\nu^2}
\begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{bmatrix},
\qquad \sigma_z = 0$$

**Estado plano de deformação (EPD):**

$$D = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix} 1-\nu & \nu & 0 \\ \nu & 1-\nu & 0 \\ 0 & 0 & \frac{1-2\nu}{2} \end{bmatrix},
\qquad \sigma_z = \nu(\sigma_x + \sigma_y)$$

Quando a saída é polar, o tensor é rotacionado para as direções radial e circunferencial pelo ângulo do centroide do elemento. Isso permite comparar diretamente $\sigma_r$ e $\sigma_\theta$ com a solução de Lamé.

---

## 3. Como compilar e rodar

### 3.1 Pré-requisitos
- `gfortran` (o código usa sintaxe Fortran legada)

### 3.2 Compilação
```bash
gfortran -std=legacy -o MEF-grupo07 MEF-grupo07.f
```

Para as malhas mais refinadas (estudo de convergência), pode ser preciso aumentar o parâmetro `npos` no início do arquivo — ele define o tamanho do vetor de trabalho. Se aparecer a mensagem de memória insuficiente, é esse valor que precisa crescer.

### 3.3 Execução

O programa pede três nomes de arquivo, em ordem: dados de entrada, saída em texto e saída em VTK.

```bash
./MEF-grupo07
```
```
Arquivo de dados: Exemplo_1.1.dat
Arquivo de saida: saida1.txt
Arquivo VTK de saida: resultado.vtk
```

### 3.4 Formato do arquivo `.dat`

A primeira linha traz os parâmetros da malha e termina com a flag `isai`, que define o sistema de saída (1 = cartesiano, 2 = polar):
```
nnode numel numat nen ndf ndm isai
```

O bloco de pressão vem depois das condições de contorno e das forças concentradas, uma linha por aresta carregada, terminando com o sentinela `0 0 0.0`:
```
ni   nj   p
```
Os nós vão na ordem anti-horária do contorno; pressão de compressão entra negativa.

---

## 4. Exemplo

`caso1_uniaxial.dat` — patch test uniaxial:

- Malha retangular 2,0 × 1,0 (T3, estado plano de tensão), 15 nós, 16 elementos.
- Apoio tipo rolete na face esquerda (nós 1, 6, 11 com $u_x=0$; o nó 1 também com $u_y=0$, apenas para eliminar movimento de corpo rígido sem restringir o efeito de Poisson).
- Pressão $p=100$ (tração) aplicada na face direita, arestas 5→10 e 10→15.

**Trecho do `.dat`:**
```
 5   10   100.00
10   15   100.00
 0    0     0.0
```

**Resultado** (idêntico nos 16 elementos):

| Campo | Esperado | Obtido |
|---|---|---|
| $\sigma_x$ | 100,0 | 100,0 |
| $\sigma_y$ | 0 | ≈ 0 |
| $\tau_{xy}$ | 0 | ≈ 0 |

Abrindo o `.vtk` no ParaView e colorindo por `sigma_x` (`CELL_DATA`), a peça aparece com cor uniforme — confirmando o campo de tensão constante.

---

## 5. Verificação

A implementação foi testada em três casos de exigência crescente:

1. **Patch test uniaxial** — pressão em uma única face. $\sigma_x = 100$ exato em todos os elementos, $\sigma_y \approx \tau_{xy} \approx 0$. Valida a conversão de pressão em força nodal equivalente.

2. **Patch test biaxial** — pressão em duas faces perpendiculares. $\sigma_x = \sigma_y = 100$ exato, $\tau_{xy} \approx 0$. Valida a soma das contribuições de arestas diferentes que compartilham um nó.

3. **Cilindro de Lamé** — quarto de anel ($a=1$, $b=2$) sob pressão interna, comparado com a solução analítica. O Q4 reproduz $\sigma_r$ e $\sigma_\theta$ com erro da ordem de 0,1%; o T3 fica na faixa esperada para um elemento de deformação constante. Refinando a malha, o erro cai de forma ordenada (ordem $\approx 2$ para o Q4), confirmando a convergência para a solução exata.

Após cada alteração do código, os resultados foram comparados com a versão anterior para garantir que nada quebrava antes de seguir em frente.

---

## 6. Referências

- Notas de Aula do Prof. Fernando L. B. Ribeiro — *Introdução ao Método dos Elementos Finitos*.
- Documentação do formato VTK legado (ASCII), para a estrutura dos blocos `POINTS`, `CELLS`, `CELL_TYPES` e `CELL_DATA`.
