# FEM-course
Códigos da disciplina de Elementos Finitos

Bem-vindo(a) ao repositório da disciplina de **Método dos Elementos Finitos**.

Aqui você encontrará os códigos e exemplos utilizados ao longo do curso.

---

## 📋 Sobre a Disciplina

Este repositório reúne os códigos e exemplos utilizados ao longo do semestre. A proposta é que você não apenas **use** os programas, mas entenda cada linha — e ao final, seja capaz de **extendê-los**.

Trabalharemos com Fortran, linguagem ainda amplamente utilizada em computação científica e engenharia. A escolha é intencional: o contato com uma linguagem compilada, eficiente e direta ajuda a entender o que acontece "por baixo dos panos" em qualquer análise FEM.

---

## 🚀 Primeiros Passos

### 1. Clone o repositório

```bash
git clone https://github.com/anabeatrizgonzaga/FEM-course.git
cd FEM-course
```

### 2. Compile o código de treliça 2D


### 3. Rode o primeiro exemplo


> **Novo no Git?** Leia o [Guia de Git para Iniciantes](guiagit.md) antes de continuar.

---

## 📐 Módulos de Código

### `01_truss2d` — Treliça 1D
Implementação clássica do método dos elementos finitos para estruturas de treliça bidimensionais. Inclui:
- Montagem da matriz de rigidez global
- Aplicação de condições de contorno
- Resolução do sistema linear
- Cálculo de deslocamento, tensões e reações de apoio

---

### `02_elastic2d` — Elasticidade Mecânica 2D
Código estruturado para análise plana (estado plano de tensão ou deformação). Inclui:
- Elemento CST (Constant Strain Triangle)
- Pós-processamento de tensões


---

## 🎯 Trabalho Final

Nas últimas semanas do semestre, cada aluno deverá **implementar uma extensão** no código `02_elastic2d`. As implementações serão definidas no decorrer do curso.

---

## 📚 Referências Principais

- **Zienkiewicz & Taylor** — *The Finite Element Method* (qualquer edição)
- **Cook et al.** — *Concepts and Applications of Finite Element Analysis*
- **Hughes** — *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*
- **Bathe** — *Finite Element Procedures*

Veja a lista completa em [`material/referencias.md`](material/referencias.md).

---

## 🤝 Como Contribuir (para os alunos)

1. Faça um **fork** do repositório
2. Crie uma branch com seu nome: `git checkout -b feat/nome-sobrenome`
3. Desenvolva sua implementação 
4. Abra um **Pull Request** descrevendo o que foi implementado

---

## 📬 Contato

**Professor:** Ana Beatriz de Carvalho Gonzaga e Silva
**E-mail:** anabeatrizgonzaga@poli.ufrj.br

---

<p align="center">
  Desenvolvido com fins didáticos · Universidade Federal do Rio de Janeiro · Departamento de Estruturas - Escola Politécnica
</p>

