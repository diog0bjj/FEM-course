c     Programa MEF - Metodo dos Elementos Finitos
c     Elasticidade Plana 2D (Estado Plano de Tensoes e Deformacoes)
c
c     Elementos disponiveis:
c       1 - T3  Estado Plano de Deformacoes (EPD)
c       2 - T3  Estado Plano de Tensoes    (EPT)
c       3 - Q4  Estado Plano de Deformacoes (EPD)
c       4 - Q4  Estado Plano de Tensoes    (EPT)
c       5 - T3  Problema de potencial
c       6 - T3  Axissimetrico 
c
c     Solver: Gradiente Conjugado Precondicionado (PCG)
c             com precondicionador diagonal (Jacobi)
c
c     Armazenamento da matriz de rigidez: Skyline (perfil de coluna)
c ======================================================================
      program mef
      parameter (npos = 500000)
      common m(npos)
      common /size/ max
      real*8 a(1)
      equivalence (m(1),a(1))
      max = npos/2
      call contr(m,a)
      stop
      end

c ======================================================================
      subroutine contr(m,a)
c     Subrotina controladora principal.
c     Orquestra todas as etapas do calculo MEF.
c ======================================================================
      common /size/ max
      integer m(*)
      real*8  a(1),dot
      character*80 fname
      external dot
c
c     Abertura de arquivos:
c
      nin  = 1
      nout = 2
      print*, 'Arquivo de dados:'
      read(*,'(a)') fname
      open(nin,file=fname)
      print*, 'Arquivo de saida:'
      read(*,'(a)') fname
      open(nout,file=fname)
c
c     Leitura das variaveis principais:
c       nnode = numero de nos
c       numel = numero de elementos
c       numat = numero de materiais
c       nen   = nos por elemento (3=triangulo, 4=quadrilatero)
c       ndf   = graus de liberdade por no (2 para 2D)
c       ndm   = dimensao do problema (2 para 2D)
c
      read(nin,*) nnode,numel,numat,nen,ndf,ndm
c
c     Verificacao basica de consistencia dos parametros de entrada:
c
      if (nnode .le. 0 .or. numel .le. 0) then
        print*, '*** ERRO: nnode e numel devem ser positivos.'
        stop
      endif
      if (nen .ne. 3 .and. nen .ne. 4) then
        print*, '*** ERRO: nen deve ser 3 (T3) ou 4 (Q4). Lido:', nen
        stop
      endif
      if (ndm .ne. 2) then
        print*, '*** ERRO: programa valido apenas para 2D (ndm=2)'
        stop
      endif
c
c     Alocacao de memoria (vetor unico m/a):
c
c      1       i1       i2       i3       i4     i5       i6
c      ---------------------------------------------------------
c      |   e   |   ie   |   ix   |   id   |   x   |   f   |
c      ---------------------------------------------------------
c
      nen1 = nen+1
      i1 = ( 1 + numat*10 - 1)*2 + 1
      i2 =  i1 + numat
      i3 =  i2 + numel*nen1
      i4 = (i3 + nnode*ndf)/2 + 1
      i5 =  i4 + nnode*ndm
      i5b = i5 + nnode*ndm
      i6 =  i5b + nnode*ndf
      call mem(i6)
c
c     Leitura dos dados:
c
      call rdata(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i5b),nnode,numel,
     .           numat,nen,ndm,ndf,nin)
c
c     Numeracao das equacoes livres (nos nao restritos):
c
      call numeq(m(i3),nnode,ndf,neq)
      print*, 'Numero de equacoes (graus de liberdade livres):', neq
c
c     Alocacao de memoria para jdiag e am:
c
c      i4    i5    i6    i7       i8
c      ------------------------------
c      |  x  |  f  |  u  | jdiag  |
c      ------------------------------
c
      i7 = (i6 + neq - 1)*2 + 1
      i8 = (i7 + neq)/2 + 1
      call mem(i8)
c
c     Perfil (skyline) da matriz de rigidez:
c
      call profil(m(i2),m(i3),m(i7),numel,nen,ndf,neq,ncs)
      print*, 'Tamanho do perfil skyline (ncs):', ncs
c
c     Alocacao de memoria para arrays locais e matriz global:
c
c      i7     i8   i9     i10    i11    i12   i13   i14
c      -----------------------------------------------
c      | jdiag | am |  xl  |  ul  |  fl  | sl  | ld  |
c      -----------------------------------------------
c
      nst = nen*ndf
      i9  =  i8  + ncs
      i10 =  i9  + nen*ndm
      i11 =  i10 + nst
      i12 =  i11 + nst
      i13 = (i12 + nst*nst - 1)*2 + 1
      i14 = (i13 + nst)/2 + 1
      call mem(i14)
c
c     Alocacao de memoria para vetores do solver PCG:
c       x=solucao, r=residuo, z=precond*r, prec=precond, d=direcao
c
c      i14   i15   i16   i17   i18   i19
c      -----------------------------------
c      |  u  |  x  |  r  |  z  | prec| d  |
c      -----------------------------------
c
      i15 = i14 + neq
      i16 = i15 + neq
      i17 = i16 + neq
      i18 = i17 + neq
      i19 = i18 + neq
c
c     Vetores de tensao por elemento (para o CELL_DATA do VTK):
c
      i20 = i19 + numel
      i21 = i20 + numel
      i22 = i21 + numel
      call mem(i22)
c
c     Forcas nodais equivalentes:
c       Transfere forcas prescritas para o vetor global u
c
      call pload(m(i3),a(i5),a(i6),nnode,ndf)
c
c     Montagem da matriz de rigidez global K e correcao do vetor f:
c       isw=2  => calcula rigidez e forcas internas do elemento
c       afl=.true. => corrige o vetor de forcas (deslocamentos prescritos)
c       bfl=.true. => monta a matriz de rigidez global
c
      call pform(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i6),m(i7),a(i8),
     .           a(i9),a(i10),a(i11),a(i12),m(i13),numel,ndm,ndf,nen,
     .           nst,2,.true.,.true.)
c
c     Resolucao do sistema de equacoes K*u = f
c     via Gradiente Conjugado Precondicionado (PCG):
c
      print*, 'Iniciando solver PCG...'
      call pcg(a(i8),a(i6),m(i7),a(i15),a(i16),a(i17),a(i18),
     .         a(i19),neq)
c
c     Impressao dos resultados (formato texto):
c
      call wdata(m(i2),m(i3),a(i4),a(i5),a(i6),nnode,numel,ndm,nen,
     .           ndf,nout)
c
c     Calculo e impressao das tensoes por elemento (T3 apenas):
c
      call tensao(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i6),a(i20),
     .            a(i21),a(i22),nnode,numel,ndm,ndf,nen,nout)
c
c     Impressao dos resultados (formato VTK para ParaView):
c
      nvtk = 3
      print*, 'Arquivo VTK de saida (ex: resultado.vtk):'
      read(*,'(a)') fname
      open(nvtk,file=fname,status='replace',action='write')
      call wvtk(m(i2),m(i3),a(i4),a(i5),a(i6),a(i20),a(i21),a(i22),
     .          nnode,numel,ndm,nen,ndf,nvtk)
      close(nvtk)
      print*, 'Arquivo VTK gerado com sucesso.'
      print*, 'Calculo concluido. Resultados gravados com sucesso.'
      return
      end

c ======================================================================
      subroutine rdata(e,ie,ix,id,x,f,feqv,nnode,numel,numat,nen,ndm,
     .                 ndf,nin)
c     Leitura dos dados de entrada:
c       e  = propriedades dos materiais  (real*8, 10 x numat)
c       ie = tipo de elemento por material (integer, numat)
c       ix = conectividade + material    (integer, nen+1 x numel)
c       id = condicoes de contorno       (integer, ndf x nnode)
c       x  = coordenadas nodais          (real*8, ndm x nnode)
c       f  = forcas nodais               (real*8, ndf x nnode)
c       feqv = forças nodais equivalentes
c ======================================================================
      integer nnode,numel,numat,ndm,nen,ndf,ie(*),ix(nen+1,*),id(ndf,*)
      integer iaux(6)
      real*8 e(10,*),x(ndm,*),f(ndf,*),feqv(ndf,*),aux(6),dum(1)
c
c     Inicializa id (todos livres) e f (sem forcas):
c
      do 5 n = 1, nnode
        do 5 j = 1, ndf
          id(j,n) = 0
          f(j,n)  = 0.d0
          feqv(j,n) = 0.d0
    5 continue
c
c     Leitura das propriedades dos materiais:
c       Linha: ma iel
c       Em seguida, a rotina do elemento le suas constantes fisicas.
c
      do 100 i = 1, numat
        read(nin,*) ma,iel
        if (ma .lt. 1 .or. ma .gt. numat) then
          print*, '*** ERRO rdata: indice de material invalido:', ma
          stop
        endif
        call elmlib(e(1,ma),dum,dum,dum,dum,1,iel,1,1,1,nin)
        ie(ma) = iel
  100 continue
c
c     Leitura das coordenadas nodais:
c       Linha: no  x  y  (z e lido mas ignorado para 2D)
c
      do 200 i = 1, nnode
        read(nin,*) k, (x(j,k),j=1,ndm)
  200 continue
c
c     Leitura da conectividade dos elementos:
c       Linha: nel  n1 n2 n3 [n4]  ma
c       Importante: nos em ordem anti-horaria!
c
      do 300 i = 1, numel
        read(nin,*) k, (ix(j,k),j=1,nen+1)
  300 continue
c
c     Leitura das condicoes de contorno (deslocamentos prescritos):
c       Linha: no  id1  id2   (id=1 => grau de liberdade restrito)
c       Termina com: 0  0  0
c
  400 continue
        read(nin,*) k, (iaux(j), j = 1, ndf)
        if (k .le. 0) goto 500
        do 410 j = 1, ndf
          id(j,k) = iaux(j)
  410   continue
      goto 400
  500 continue
c
c     Leitura das forcas nodais aplicadas:
c       Linha: no  fx  fy
c       Termina com: 0  0.0  0.0
c
  510 continue
        read(nin,*) k, (aux(j), j = 1, ndf)
        if (k .le. 0) goto 600
        do 520 j = 1, ndf
          f(j,k) = aux(j)
  520   continue
      goto 510
  600 continue
c
c     Calcula forcas equivalentes de pressao e soma em f:
c
      call pressao(x,feqv,ix,ndm,ndf,nnode,nen,numel,nin)
      do 610 n = 1, nnode
        do 610 j = 1, ndf
          if (id(j,n) .eq. 0) f(j,n) = f(j,n) + feqv(j,n)
  610 continue
      return
      end

c ======================================================================
      subroutine pressao(x,feqv,ix,ndm,ndf,nnode,nen,numel,nin)
c     Leitura da pressao distribuida no contorno e conversao em
c     forcas nodais equivalentes.
c       Linha de entrada: no_i  no_j  pressao
c       Termina com: 0  0  0.0
c
c     Convencao de sinal: pressao positiva = tracao (normal aponta
c     para fora do dominio). Os nos i,j devem ser fornecidos na
c     mesma ordem de percurso do contorno externo (anti-horario),
c     ou seja: sentido em que o solido fica sempre a esquerda de
c     quem caminha de ni para nj.
c
c     Cada aresta (ni,nj) e' automaticamente localizada dentro da
c     malha: procura-se, em ix, o elemento (T3 ou Q4) que possui
c     ni e nj como nos consecutivos de sua conectividade local.
c     Isso serve tanto para validar a entrada (acusa erro se a
c     aresta nao existir) quanto para confirmar que o sentido dado
c     e' consistente com a orientacao anti-horaria do elemento.
c ======================================================================
      integer ndm,ndf,nnode,nen,numel,nin,ni,nj
      integer ix(nen+1,*)
      real*8  x(ndm,*),feqv(ndf,*)
      real*8  p
      integer nelem,isense
c
  10  continue
        read(nin,*) ni, nj, p
        if (ni .le. 0) goto 100
c
c       Identificacao automatica do elemento/aresta (ni,nj):
c
        call findedge(ix,nen,numel,ni,nj,nelem,isense)
        if (nelem .eq. 0) then
          print*, '*** ERRO pressao: aresta (',ni,',',nj,
     .            ') nao pertence a nenhum elemento da malha.'
          print*, '    Verifique os numeros de no informados.'
          stop
        endif
        if (isense .lt. 0) then
          print*, '*** AVISO pressao: aresta (',ni,',',nj,
     .            ') do elemento',nelem,
     .            ' foi lida no sentido invertido em relacao'
          print*, '    a conectividade anti-horaria do elemento.'
          print*, '    O sinal da carga equivalente sera invertido',
     .            ' (verifique se e o pretendido).'
        endif
c
c       Forcas nodais equivalentes por integracao consistente:
c
        call edgfor(x,feqv,ndm,ndf,ni,nj,p)
      goto 10
 100  continue
      return
      end

c ======================================================================
      subroutine findedge(ix,nen,numel,ni,nj,nelem,isense)
c     Localiza, entre todos os elementos, aquele cuja conectividade
c     local contem a aresta (ni,nj) como par de nos consecutivos
c     (identificacao automatica da aresta pertencente ao elemento).
c
c     Saida:
c       nelem  = numero do elemento encontrado (0 = nao encontrado)
c       isense = +1  se (ni,nj) aparece na mesma ordem da
c                    conectividade local (anti-horaria, como lida)
c                -1  se aparece invertida, (nj,ni)
c
c     Percorre as nen arestas locais de cada elemento:
c       T3 (nen=3): (v1,v2),(v2,v3),(v3,v1)
c       Q4 (nen=4): (v1,v2),(v2,v3),(v3,v4),(v4,v1)
c ======================================================================
      integer nen,numel,ni,nj,nelem,isense
      integer ix(nen+1,*)
      integer nel,k,k2,v1,v2
c
      nelem  = 0
      isense = 0
      do 200 nel = 1, numel
        do 100 k = 1, nen
          k2 = mod(k,nen) + 1
          v1 = ix(k ,nel)
          v2 = ix(k2,nel)
          if (v1 .eq. ni .and. v2 .eq. nj) then
            nelem  = nel
            isense = 1
            return
          else if (v1 .eq. nj .and. v2 .eq. ni) then
            nelem  = nel
            isense = -1
            return
          endif
  100   continue
  200 continue
      return
      end

c ======================================================================
      subroutine edgfor(x,feqv,ndm,ndf,ni,nj,p)
c     Forcas nodais equivalentes de uma pressao p, uniforme ao longo
c     da aresta linear (ni,nj), obtidas por integracao consistente
c     com as MESMAS funcoes de forma lineares usadas pelo elemento
c     (T3 e Q4 tem aresta de 2 nos, interpolacao linear em ambos):
c
c       N1(csi) = (1-csi)/2 ,  N2(csi) = (1+csi)/2 ,  csi em [-1,1]
c
c     Vetor de tracao (traction) na aresta, com a normal saindo do
c     dominio (convencao: solido a esquerda de ni->nj):
c       t = p * n ,   n = (dy,-dx)/L  (dx,dy = componentes de nj-ni)
c
c     Forca nodal consistente:
c       F_a = integral( N_a(csi) * t ) * (L/2) dcsi ,  a = 1,2
c     calculada aqui por quadratura de Gauss de 2 pontos (exata para
c     qualquer integrando ate grau 3; aqui N_a*t e' grau 1, portanto
c     o resultado e EXATO, igual ao da formula fechada F=(p L/2)*n,
c     mas agora obtido explicitamente a partir de N1,N2 - preparado
c     para o caso de p variar ao longo da aresta (ex.: pressao
c     hidrostatica) bastando interpolar p(csi)=N1*p1+N2*p2.
c ======================================================================
      integer ndm,ndf,ni,nj
      real*8  x(ndm,*),feqv(ndf,*),p
      real*8  dx,dy,L,nx,ny,tx,ty,detj
      real*8  gp(2),gw(2),csi,n1,n2,f1,f2
      integer ig
c
      dx = x(1,nj) - x(1,ni)
      dy = x(2,nj) - x(2,ni)
      L  = dsqrt(dx*dx + dy*dy)
      nx =  dy / L
      ny = -dx / L
      tx = p * nx
      ty = p * ny
c
c     Jacobiano da transformacao csi (em [-1,1]) -> comprimento de
c     arco s (em [0,L]): ds/dcsi = L/2 (constante, aresta reta)
c
      detj = L / 2.d0
c
c     Pontos e pesos de Gauss (2 pontos, exato ate grau 3):
c
      gp(1) = -1.d0/dsqrt(3.d0)
      gp(2) =  1.d0/dsqrt(3.d0)
      gw(1) = 1.d0
      gw(2) = 1.d0
c
      f1 = 0.d0
      f2 = 0.d0
      do 10 ig = 1, 2
        csi = gp(ig)
        n1  = (1.d0-csi)/2.d0
        n2  = (1.d0+csi)/2.d0
        f1  = f1 + gw(ig) * n1 * detj
        f2  = f2 + gw(ig) * n2 * detj
   10 continue
c
c     f1 e f2 valem L/2 cada (verificacao: integral de N1 e N2 ao
c     longo da aresta = L/2), reproduzindo a formula classica.
c
      feqv(1,ni) = feqv(1,ni) + tx * f1
      feqv(2,ni) = feqv(2,ni) + ty * f1
      feqv(1,nj) = feqv(1,nj) + tx * f2
      feqv(2,nj) = feqv(2,nj) + ty * f2
      return
      end

c ======================================================================
      subroutine tensao(e,ie,ix,id,x,f,u,sxa,sya,txya,nnode,numel,ndm,
     .                  ndf,nen,nout)
c     Calculo das tensoes em cada elemento, a partir dos deslocamentos
c     resultantes da solucao do sistema.
c       iel=1 : T3-EPD  (tensao constante no elemento)
c       iel=2 : T3-EPT  (tensao constante no elemento)
c       iel=3 : Q4-EPD  (tensao avaliada no centro do elemento,
c                        ponto de Gauss de 1 ponto, r=s=0)
c       iel=4 : Q4-EPT  (idem)
c
c     sigma = D * B * u_local   (apostila, eq. 2.41/2.6.1)
c
c     sxa, sya, txya: vetores de saida (tamanho numel) com a tensao
c     de cada elemento, para uso posterior na escrita do VTK (CELL_DATA)
c ======================================================================
      integer nnode,numel,ndm,ndf,nen,nout
      integer ie(*),ix(nen+1,*),id(ndf,*)
      real*8  e(10,*),x(ndm,*),f(ndf,*),u(*)
      real*8  sxa(*),sya(*),txya(*)
      integer nel,ma,iel,i,no,k,nenel
      real*8  x1,y1,x2,y2,x3,y3,det
      real*8  xji11,xji12,xji21,xji22
      real*8  hx(4),hy(4),hr(4),hs(4),ul(2,4)
      real*8  xj11,xj12,xj21,xj22
      real*8  epsx,epsy,gamxy
      real*8  my,nu,a,b,d11,d12,d22,d33
      real*8  sx,sy,txy
c
      write(nout,'(a)') 'Tensoes por elemento (nel, sx, sy, txy):'
      do 500 nel = 1, numel
        ma  = ix(nen+1,nel)
        iel = ie(ma)
c
c       So calcula para T3 ou Q4 (iel=1,2,3,4):
c
        if (iel .lt. 1 .or. iel .gt. 4) then
          sxa(nel)  = 0.d0
          sya(nel)  = 0.d0
          txya(nel) = 0.d0
          goto 500
        endif
        if (iel .le. 2) then
          nenel = 3
        else
          nenel = 4
        endif
c
c       Deslocamento local de cada no do elemento, buscado da
c       solucao global (mesma logica que wdata ja usa: se o grau
c       de liberdade eh restrito, id=0 e o valor vem de f(); se eh
c       livre, id>0 aponta pra posicao na solucao u()):
c
        do 100 i = 1, nenel
          no = ix(i,nel)
          do 100 k = 1, ndf
            ul(k,i) = f(k,no)
            if (id(k,no) .gt. 0) ul(k,i) = u(id(k,no))
  100   continue
c
        if (iel .le. 2) then
c
c         --- T3: deformacao constante (identico a elmt01/elmt02) ---
c
          x1 = x(1,ix(1,nel))
          y1 = x(2,ix(1,nel))
          x2 = x(1,ix(2,nel))
          y2 = x(2,ix(2,nel))
          x3 = x(1,ix(3,nel))
          y3 = x(2,ix(3,nel))
c
          det   = (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
          xji11 =  (y2-y3)/det
          xji12 = -(y1-y3)/det
          xji21 = -(x2-x3)/det
          xji22 =  (x1-x3)/det
          hx(1) =  xji11
          hx(2) =  xji12
          hx(3) = -hx(1)-hx(2)
          hy(1) =  xji21
          hy(2) =  xji22
          hy(3) = -hy(1)-hy(2)
c
          epsx  = hx(1)*ul(1,1) + hx(2)*ul(1,2) + hx(3)*ul(1,3)
          epsy  = hy(1)*ul(2,1) + hy(2)*ul(2,2) + hy(3)*ul(2,3)
          gamxy = hy(1)*ul(1,1) + hx(1)*ul(2,1)
     .          + hy(2)*ul(1,2) + hx(2)*ul(2,2)
     .          + hy(3)*ul(1,3) + hx(3)*ul(2,3)
c
        else
c
c         --- Q4: avaliado no centro do elemento (r=s=0), que
c             equivale ao ponto de integracao de Gauss de 1 ponto
c             (apostila, Tabela 7.1: n=1, xi=0). Mesmas formulas
c             de elmt03/elmt04, so que particularizadas em r=s=0 ---
c
          hr(1) =  0.25d0
          hr(2) = -0.25d0
          hr(3) = -0.25d0
          hr(4) =  0.25d0
          hs(1) =  0.25d0
          hs(2) =  0.25d0
          hs(3) = -0.25d0
          hs(4) = -0.25d0
c
          xj11 = hr(1)*x(1,ix(1,nel))+hr(2)*x(1,ix(2,nel))
     .         + hr(3)*x(1,ix(3,nel))+hr(4)*x(1,ix(4,nel))
          xj21 = hs(1)*x(1,ix(1,nel))+hs(2)*x(1,ix(2,nel))
     .         + hs(3)*x(1,ix(3,nel))+hs(4)*x(1,ix(4,nel))
          xj12 = hr(1)*x(2,ix(1,nel))+hr(2)*x(2,ix(2,nel))
     .         + hr(3)*x(2,ix(3,nel))+hr(4)*x(2,ix(4,nel))
          xj22 = hs(1)*x(2,ix(1,nel))+hs(2)*x(2,ix(2,nel))
     .         + hs(3)*x(2,ix(3,nel))+hs(4)*x(2,ix(4,nel))
c
          det   = xj11*xj22 - xj21*xj12
          xji11 =  xj22/det
          xji12 = -xj12/det
          xji21 = -xj21/det
          xji22 =  xj11/det
c
          do 200 k = 1, 4
            hx(k) = xji11*hr(k) + xji12*hs(k)
            hy(k) = xji21*hr(k) + xji22*hs(k)
  200     continue
c
          epsx  = hx(1)*ul(1,1)+hx(2)*ul(1,2)+hx(3)*ul(1,3)
     .          + hx(4)*ul(1,4)
          epsy  = hy(1)*ul(2,1)+hy(2)*ul(2,2)+hy(3)*ul(2,3)
     .          + hy(4)*ul(2,4)
          gamxy = hy(1)*ul(1,1)+hx(1)*ul(2,1)
     .          + hy(2)*ul(1,2)+hx(2)*ul(2,2)
     .          + hy(3)*ul(1,3)+hx(3)*ul(2,3)
     .          + hy(4)*ul(1,4)+hx(4)*ul(2,4)
        endif
c
c       Matriz constitutiva D (apostila, eq. 2.8 EPT / eq. 2.9 EPD):
c       EPD: iel=1 (T3) ou iel=3 (Q4). EPT: iel=2 (T3) ou iel=4 (Q4).
c
        my = e(1,ma)
        nu = e(2,ma)
        if (iel .eq. 1 .or. iel .eq. 3) then
c         Estado Plano de Deformacao (EPD):
          a   = 1.d0+nu
          b   = a*(1.d0-2.d0*nu)
          d11 = my*(1.d0-nu)/b
          d12 = my*nu/b
          d22 = d11
          d33 = my/(2.d0*a)
        else
c         Estado Plano de Tensao (EPT):
          d11 = my/(1.d0-nu*nu)
          d12 = d11*nu
          d22 = d11
          d33 = my/(2.d0*(1.d0+nu))
        endif
c
c       Tensoes: sigma = D * eps
c
        sx  = d11*epsx + d12*epsy
        sy  = d12*epsx + d22*epsy
        txy = d33*gamxy
c
        sxa(nel)  = sx
        sya(nel)  = sy
        txya(nel) = txy
c
        write(nout,'(i8,3e15.5)') nel, sx, sy, txy
  500 continue
      return
      end

c ======================================================================
      subroutine numeq(id,nnode,ndf,neq)
c     Numeracao sequencial das equacoes livres.
c     Nos livres  (id=0) recebem numero de equacao: 1, 2, 3, ...
c     Nos restritos (id=1) recebem id=0 (sem equacao associada).
c ======================================================================
      integer nnode,ndf,neq,id(ndf,*)
      neq = 0
      do 100 n = 1, nnode
      do 100 i = 1, ndf
        j = id(i,n)
        if (j .eq. 0) then
          neq = neq + 1
          id(i,n) = neq
        else
          id(i,n) = 0
        endif
  100 continue
      return
      end

c ======================================================================
      subroutine profil(ix,id,jdiag,numel,nen,ndf,neq,ncs)
c     Calculo do perfil (skyline) da matriz de rigidez.
c     Para cada coluna j, determina a altura efetiva (numero de
c     coeficientes nao-nulos acima da diagonal), com base na
c     conectividade dos elementos.
c     Preenche jdiag(i) = posicao da diagonal do coeficiente i
c     no vetor de armazenamento skyline.
c ======================================================================
      integer numel,nen,ndf,neq,ncs,ix(nen+1,*),id(ndf,*),jdiag(*)
      call mzero(jdiag,1,neq)
c
c     Alturas de coluna: para cada par de equacoes (kk, ll)
c     pertencentes ao mesmo elemento, a altura da coluna
c     max(kk,ll) deve incluir o coeficiente |kk-ll|.
c
      do 400 nel = 1, numel
      do 400 i = 1, nen
        noi = ix(i,nel)
        if (noi .eq. 0) goto 400
        do 300 k = 1, ndf
          kk = id(k,noi)
          if (kk .eq. 0) goto 300
          do 200 j = i, nen
            noj = ix(j,nel)
            if (noj .eq. 0) goto 200
            do 100 l = 1, ndf
              ll = id(l,noj)
              if (ll .eq. 0) goto 100
              m = max0(kk,ll)
              jdiag(m) = max0(jdiag(m),iabs(kk-ll))
  100       continue
  200     continue
  300   continue
  400 continue
c
c     Ponteiros da diagonal: converte alturas em posicoes absolutas.
c     jdiag(i) = jdiag(i-1) + altura(i) + 1
c
      ncs = 1
      jdiag(1) = 1
      if (neq .eq. 1) return
      do 500 i = 2, neq
        jdiag(i) = jdiag(i) + jdiag(i-1) + 1
  500 continue
      ncs = jdiag(neq)
      return
      end

c ======================================================================
      subroutine pload(id,f,u,nnode,ndf)
c     Montagem do vetor de forcas nodais globais.
c     Soma as equivalentes
c     Transfere as forcas do array local f(ndf,nnode) para o
c     vetor global u(neq), usando a numeracao de equacoes id.
c ======================================================================
      integer nnode, ndf, id(ndf,*)
      real*8 f(ndf,*), u(*)
      do 100 i = 1, nnode
      do 100 j = 1, ndf
        k = id(j,i)
        if (k .gt. 0) u(k) = f(j,i)
  100 continue
      return
      end

c ======================================================================
      subroutine pform(e,ie,ix,id,x,f,u,jdiag,am,xl,ul,fl,sl,ld,
     .                 numel,ndm,ndf,nen,nst,isw,afl,bfl)
c     Loop nos elementos: monta a matriz de rigidez global e
c     corrige o vetor de forcas para deslocamentos prescritos.
c       isw = codigo de instrucao para a rotina de elemento
c       afl = .true. => corrige o vetor de forcas
c       bfl = .true. => monta a matriz de rigidez global
c ======================================================================
      integer numel,ndm,ndf,nen,nst,isw
      integer ie(*),ix(nen+1,*),id(ndf,*),jdiag(*),ld(ndf,*)
      real*8  e(10,*),x(ndm,*),f(ndf,*),u(*),am(*),xl(ndm,*),ul(ndf,*)
      real*8  fl(*),sl(*)
      logical afl, bfl
c
c     Loop nos elementos:
c
      do 700 nel = 1, numel
c
c       Loop nos nos do elemento: monta vetores locais xl, ul, ld
c
        do 600 i = 1, nen
          no = ix(i,nel)
          if (no .gt. 0) goto 300
c
c         No = 0: zera vetores locais
c
          do 100 j = 1, ndm
            xl(j,i) = 0.d0
  100     continue
          do 200 j = 1, ndf
            ul(j,i) = 0.d0
            ld(j,i) = 0
  200     continue
          goto 600
c
c         No > 0: copia coordenadas e deslocamentos/forcas prescritas
c
  300     continue
          do 400 j = 1, ndm
            xl(j,i) = x(j,no)
  400     continue
          do 500 j = 1, ndf
c           Numeracao global da equacao do grau de liberdade (j,no):
            k = id(j,no)
            ld(j,i) = k
c           Deslocamento prescrito (k=0 indica no restrito):
            ul(j,i) = 0.d0
            if (k .eq. 0) ul(j,i) = f(j,no)
  500     continue
  600   continue
        ma  = ix(nen+1,nel)
        iel = ie(ma)
c
c       Calcula rigidez e forcas locais do elemento:
c
        call elmlib(e(1,ma),xl,ul,fl,sl,nel,iel,ndm,nst,isw,1)
c
c       Acumula contribuicoes local -> global:
c
        call addstf(am,u,jdiag,fl,sl,ld,nst,afl,bfl)
  700 continue
      return
      end

c ======================================================================
      subroutine addstf(a,b,jdiag,p,s,ld,nst,afl,bfl)
c     Montagem local->global da matriz de rigidez e vetor de forcas.
c     Armazenamento skyline: apenas a parte triangular inferior
c     de cada coluna e armazenada.
c       p = vetor de forcas locais do elemento
c       s = matriz de rigidez local do elemento
c       ld = numeracao global das equacoes do elemento
c ======================================================================
      integer jdiag(*),ld(*),nst
      real*8 a(*),b(*),p(*),s(nst,*)
      logical afl,bfl
      do 200 j = 1, nst
        k = ld(j)
        if (k .eq. 0) goto 200
c
c       Correcao do vetor de forcas para deslocamentos prescritos:
c       b(k) -= s(i,j) * u_prescrito  (contribuicao de cada gdl)
c
        if (afl) b(k) = b(k) - p(j)
        if (.not. bfl) goto 200
c
c       Montagem na matriz de rigidez (parte triangular inferior):
c       l = offset base da coluna k no vetor skyline
c
        l = jdiag(k) - k
        do 100 i = 1, nst
          m = ld(i)
          if (m .gt. k .or. m .eq. 0) goto 100
          m = l + m
          a(m) = a(m) + s(i,j)
  100   continue
  200 continue
      return
      end

c ======================================================================
      subroutine elmlib(e,xl,ul,fl,sl,nel,iel,ndm,nst,isw,nin)
c     Biblioteca de elementos: despachante para cada tipo de elemento.
c       iel = 1 => T3 Estado Plano de Deformacoes (EPD)
c       iel = 2 => T3 Estado Plano de Tensoes    (EPT)
c       iel = 3 => Q4 Estado Plano de Deformacoes (EPD)
c       iel = 4 => Q4 Estado Plano de Tensoes    (EPT)
c       iel = 5 => T3 Conducao de calor 2D
c       iel = 6 => T3 Axissimetrico
c ======================================================================
      integer nel,iel,ndm,nst,isw
      real*8 e(*),xl(*),ul(*),fl(*),sl(*)
      goto (100,200,300,400,500,600) iel
      print*, '*** ERRO elmlib: tipo de elemento ',iel,
     .        ' nao existe. Elemento: ',nel
      stop
  100 call elmt01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  200 call elmt02(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  300 call elmt03(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  400 call elmt04(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  500 call elmt01_t(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
  600 call elmt01_axi(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
      end

c ======================================================================
      subroutine elmt01(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento triangular linear T3
c     Hipotese: Estado Plano de DEFORMACOES (EPD)
c     Parametros do material: E (modulo de Young), nu (Poisson)
c
c     Matriz constitutiva EPD:
c       c = E(1-nu)/[(1+nu)(1-2nu)]
c       D = | c          c*nu/(1-nu)  0            |
c           | c*nu/(1-nu)  c          0            |
c           | 0            0    E/[2(1+nu)]        |
c ======================================================================
      integer nel,ndm,nst,isw
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  my,nu,a,b,c,wt
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas do material
c
  100 continue
      read(nin,*) e(1),e(2)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz Jacobiana:
c       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
c           | x2-x3  y2-y3 |
c
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
c
c     Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao em relacao a x e y:
c       hx(i) = dN_i/dx,  hy(i) = dN_i/dy
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva D (Estado Plano de Deformacoes):
c
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
c
c     Matriz de rigidez: Ke = B^T D B * det/2
c       wt = Area do triangulo = det/2
c
      wt = 0.5d0*det
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
c         s(l,k)     += (dNi/dx*d11*dNj/dx + dNi/dy*d33*dNj/dy)*A
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
c         s(l,k+1)   += (dNi/dx*d12*dNj/dy + dNi/dy*d33*dNj/dx)*A
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
c         s(l+1,k)   += (dNi/dy*d21*dNj/dx + dNi/dx*d33*dNj/dy)*A
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
c         s(l+1,k+1) += (dNi/dy*d22*dNj/dy + dNi/dx*d33*dNj/dx)*A
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
c
c     Vetor de forcas internas: p = Ke * u
c
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end

c ======================================================================
      subroutine elmt02(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento triangular linear T3
c     Hipotese: Estado Plano de TENSOES (EPT)
c     Parametros do material: E, nu, espessura t
c
c     Matriz constitutiva EPT:
c       a = E/(1-nu^2)
c       D = | a      a*nu      0          |
c           | a*nu   a         0          |
c           | 0      0    E/[2(1+nu)]     |
c ======================================================================
      integer nel,ndm,nst,isw
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  my,nu,thic,a,wt
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas: E, nu, espessura
c
  100 continue
      read(nin,*) e(1),e(2),e(3)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz Jacobiana:
c
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
c
c     Inversa da Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva D (Estado Plano de Tensoes):
c
      my   = e(1)
      nu   = e(2)
      thic = e(3)
      a    = my/(1.d0-nu*nu)
      d11  = a
      d12  = a*nu
      d21  = d12
      d22  = a
      d33  = my/(2.d0*(1.d0+nu))
c
c     Matriz de rigidez: Ke = B^T D B * (det/2) * t
c       wt = Area * espessura
c
      wt = 0.5d0*det*thic
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt02: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end

c ======================================================================
      subroutine elmt03(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento quadrilatero bilinear Q4
c     Hipotese: Estado Plano de DEFORMACOES (EPD)
c     Parametros do material: E, nu
c
c     Integracao numerica: Gauss 2x2 (4 pontos)
c     Pontos de Gauss: (r,s) = (+/-1/sqrt(3), +/-1/sqrt(3))
c     Pesos: wi = wj = 1.0 (todos iguais, ja incorporados)
c
c     Funcoes de interpolacao bilineares (coordenadas isoparametricas):
c       N1=(1+r)(1+s)/4, N2=(1-r)(1+s)/4
c       N3=(1-r)(1-s)/4, N4=(1+r)(1-s)/4
c ======================================================================
      integer nel,ndm,nst,isw,i,j,k,l,m,n
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
      real*8  det,hx(4),hy(4),xji(ndm,2),hr(4),hs(4)
      real*8  rr,ss,d11,d12,d21,d22,d33
      real*8  my,nu,a,b,c,wt
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas: E, nu
c
  100 continue
      read(nin,*) e(1),e(2)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz constitutiva D (Estado Plano de Deformacoes):
c
  200 continue
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
c
c     Zeragem da matriz de rigidez local:
c
      do 20 i = 1, nst
      do 20 j = 1, nst
        s(i,j) = 0.d0
   20 continue
c
c     Loop nos 4 pontos de Gauss (2x2):
c     Coordenadas: +/-1/sqrt(3) = +/-0.577350269189626
c     Pesos: todos 1.0
c
      do 600 i = 1, 2
      do 600 j = 1, 2
        if (i .eq. 1 .and. j .eq. 1) then
          rr =  0.577350269189626d0
          ss = -0.577350269189626d0
        else if (i .eq. 1 .and. j .eq. 2) then
          rr =  0.577350269189626d0
          ss =  0.577350269189626d0
        else if (i .eq. 2 .and. j .eq. 1) then
          rr = -0.577350269189626d0
          ss =  0.577350269189626d0
        else
          rr = -0.577350269189626d0
          ss = -0.577350269189626d0
        endif
c
c       Derivadas das funcoes de interpolacao em r e s:
c         hr(i) = dNi/dr,  hs(i) = dNi/ds
c
        hr(1) =   (1.d0+ss) / 4.d0
        hr(2) = - (1.d0+ss) / 4.d0
        hr(3) = - (1.d0-ss) / 4.d0
        hr(4) =   (1.d0-ss) / 4.d0
        hs(1) =   (1.d0+rr) / 4.d0
        hs(2) =   (1.d0-rr) / 4.d0
        hs(3) = - (1.d0-rr) / 4.d0
        hs(4) = - (1.d0+rr) / 4.d0
c
c       Matriz Jacobiana: J = dX/d(r,s)
c         J(1,1)=sum(hr*x), J(2,1)=sum(hs*x)
c         J(1,2)=sum(hr*y), J(2,2)=sum(hs*y)
c
        xj(1,1) = hr(1)*x(1,1)+hr(2)*x(1,2)+hr(3)*x(1,3)+hr(4)*x(1,4)
        xj(2,1) = hs(1)*x(1,1)+hs(2)*x(1,2)+hs(3)*x(1,3)+hs(4)*x(1,4)
        xj(1,2) = hr(1)*x(2,1)+hr(2)*x(2,2)+hr(3)*x(2,3)+hr(4)*x(2,4)
        xj(2,2) = hs(1)*x(2,1)+hs(2)*x(2,2)+hs(3)*x(2,3)+hs(4)*x(2,4)
c
c       Determinante da Jacobiana: det(J) = J11*J22 - J21*J12
c       det > 0 garante elemento valido (numeracao anti-horaria)
c
        det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
        if (det .le. 0.d0) goto 1000
c
c       Inversa da Jacobiana: J^-1 = (1/det)*adj(J)
c
        xji(1,1) =  xj(2,2) / det
        xji(1,2) = -xj(1,2) / det
        xji(2,1) = -xj(2,1) / det
        xji(2,2) =  xj(1,1) / det
c
c       Derivadas das funcoes de interpolacao em x e y:
c         hx(k) = dNk/dx = J^-1(1,1)*hr(k) + J^-1(1,2)*hs(k)
c         hy(k) = dNk/dy = J^-1(2,1)*hr(k) + J^-1(2,2)*hs(k)
c
        do 300 k = 1, 4
          hx(k) = xji(1,1)*hr(k) + xji(1,2)*hs(k)
          hy(k) = xji(2,1)*hr(k) + xji(2,2)*hs(k)
  300   continue
c
c       Contribuicao ao Ke: Ke += B^T*D*B * det(J) * wi*wj  (wi=wj=1)
c
        wt = det
        do 500 m = 1, 4
          k = (m-1)*2+1
          do 400 n = 1, 4
            l = (n-1)*2+1
            s(l,k)     = s(l,k)    +(hx(n)*d11*hx(m)+hy(n)*d33*hy(m))*wt
            s(l,k+1)   = s(l,k+1)  +(hx(n)*d12*hy(m)+hy(n)*d33*hx(m))*wt
            s(l+1,k)   = s(l+1,k)  +(hy(n)*d12*hx(m)+hx(n)*d33*hy(m))*wt
            s(l+1,k+1) = s(l+1,k+1)+(hy(n)*d22*hy(m)+hx(n)*d33*hx(m))*wt
  400     continue
  500   continue
  600 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt03: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end

c ======================================================================
      subroutine elmt04(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento quadrilatero bilinear Q4
c     Hipotese: Estado Plano de TENSOES (EPT)
c     Parametros do material: E, nu, espessura t
c
c     Integracao numerica: Gauss 2x2 (4 pontos)
c     Igual ao elmt03, com matriz D de EPT e fator espessura.
c ======================================================================
      integer nel,ndm,nst,isw,i,j,k,l,m,n
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
      real*8  det,hx(4),hy(4),xji(ndm,2),hr(4),hs(4)
      real*8  rr,ss,d11,d12,d21,d22,d33
      real*8  my,nu,a,thic,wt
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas: E, nu, espessura
c
  100 continue
      read(nin,*) e(1),e(2),e(3)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz constitutiva D (Estado Plano de Tensoes):
c
  200 continue
      my   = e(1)
      nu   = e(2)
      thic = e(3)
      a    = my/(1.d0-nu*nu)
      d11  = a
      d12  = a*nu
      d21  = d12
      d22  = a
      d33  = my/(2.d0*(1.d0+nu))
c
c     Zeragem da matriz de rigidez local:
c
      do 20 i = 1, nst
      do 20 j = 1, nst
        s(i,j) = 0.d0
   20 continue
c
c     Loop nos 4 pontos de Gauss (2x2):
c
      do 600 i = 1, 2
      do 600 j = 1, 2
        if (i .eq. 1 .and. j .eq. 1) then
          rr =  0.577350269189626d0
          ss = -0.577350269189626d0
        else if (i .eq. 1 .and. j .eq. 2) then
          rr =  0.577350269189626d0
          ss =  0.577350269189626d0
        else if (i .eq. 2 .and. j .eq. 1) then
          rr = -0.577350269189626d0
          ss =  0.577350269189626d0
        else
          rr = -0.577350269189626d0
          ss = -0.577350269189626d0
        endif
c
c       Derivadas das funcoes de interpolacao em r e s:
c
        hr(1) =   (1.d0+ss) / 4.d0
        hr(2) = - (1.d0+ss) / 4.d0
        hr(3) = - (1.d0-ss) / 4.d0
        hr(4) =   (1.d0-ss) / 4.d0
        hs(1) =   (1.d0+rr) / 4.d0
        hs(2) =   (1.d0-rr) / 4.d0
        hs(3) = - (1.d0-rr) / 4.d0
        hs(4) = - (1.d0+rr) / 4.d0
c
c       Matriz Jacobiana:
c
        xj(1,1) = hr(1)*x(1,1)+hr(2)*x(1,2)+hr(3)*x(1,3)+hr(4)*x(1,4)
        xj(2,1) = hs(1)*x(1,1)+hs(2)*x(1,2)+hs(3)*x(1,3)+hs(4)*x(1,4)
        xj(1,2) = hr(1)*x(2,1)+hr(2)*x(2,2)+hr(3)*x(2,3)+hr(4)*x(2,4)
        xj(2,2) = hs(1)*x(2,1)+hs(2)*x(2,2)+hs(3)*x(2,3)+hs(4)*x(2,4)
c
c       Determinante da Jacobiana:
c
        det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
        if (det .le. 0.d0) goto 1000
c
c       Inversa da Jacobiana:
c
        xji(1,1) =  xj(2,2) / det
        xji(1,2) = -xj(1,2) / det
        xji(2,1) = -xj(2,1) / det
        xji(2,2) =  xj(1,1) / det
c
c       Derivadas das funcoes de interpolacao em x e y:
c
        do 300 k = 1, 4
          hx(k) = xji(1,1)*hr(k) + xji(1,2)*hs(k)
          hy(k) = xji(2,1)*hr(k) + xji(2,2)*hs(k)
  300   continue
c
c       Contribuicao ao Ke: Ke += B^T*D*B * det(J) * t * wi*wj
c
        wt = det * thic
        do 500 m = 1, 4
          k = (m-1)*2+1
          do 400 n = 1, 4
            l = (n-1)*2+1
            s(l,k)     = s(l,k)    +(hx(n)*d11*hx(m)+hy(n)*d33*hy(m))*wt
            s(l,k+1)   = s(l,k+1)  +(hx(n)*d12*hy(m)+hy(n)*d33*hx(m))*wt
            s(l+1,k)   = s(l+1,k)  +(hy(n)*d12*hx(m)+hx(n)*d33*hy(m))*wt
            s(l+1,k+1) = s(l+1,k+1)+(hy(n)*d22*hy(m)+hx(n)*d33*hx(m))*wt
  400     continue
  500   continue
  600 continue
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt04: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end

c ======================================================================
      subroutine elmt01_t(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento triangular linear T3
c     Hipotese: Conducao de calor 2D (equacao do calor) 
c     Parametros do material: k ( condutividade termica)
c
c ======================================================================
      integer nel,ndm,nst,isw
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
      real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  difus,wt
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas do material
c
  100 continue
      read(nin,*) e(3)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz Jacobiana:
c       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
c           | x2-x3  y2-y3 |
c
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      if (det .le. 0.d0) goto 1000
c
c     Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao em relacao a x e y:
c       hx(i) = dN_i/dx,  hy(i) = dN_i/dy
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
      difus  = e(3)
c
c       wt = Area do triangulo = det/2
c
      wt = 0.5d0*det
      do 220 j = 1, 3
        do 210 i = 1, 3
c         s(i,j)     += difus*(dNi/dx*dNj/dx + dNi/dy*dNj/dy)*A
          s(i,j)     = ( hx(i)*difus*hx(j) + hy(i)*difus*hy(j) ) * wt
  210   continue
  220 continue
c
c     Vetor de forcas internas: p = Ke * u
c
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01_t: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end
c ======================================================================
      subroutine elmt01_axi(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     Elemento triangular linear T3
c     Hipotese: Axissimetrico
c     Parametros do material: E (modulo de Young), nu (Poisson)
c
c     Matriz constitutiva EPD:
c       c = E(1-nu)/[(1+nu)(1-2nu)]
c       D = | c          c*nu/(1-nu)  0           c*nu/(1-nu) |     |
c           | c*nu/(1-nu)  c          0           c*nu/(1-nu) |
c           | 0            0    E/[2(1+nu)]             0     |
c           | c*nu/(1-nu)  c*nu/(1-nu)  0               c     |      
c ======================================================================
      integer nel,ndm,nst,isw
      real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
      real*8  det,xj11,xj12,xj21,xj22,hx(4),hy(4)
      real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
      real*8  my,nu,a,b,c,wt
      real*8   PI
      parameter (PI = 3.1415926535)
      goto (100,200) isw
c
c     isw=1: leitura das constantes fisicas do material
c
  100 continue
      read(nin,*) e(1),e(2)
      return
c
c     isw=2: calculo da matriz de rigidez
c
c     Matriz Jacobiana:
c       J = | x1-x3  y1-y3 |    det(J) = 2 * Area do triangulo
c           | x2-x3  y2-y3 |
c
  200 continue
      xj11 = x(1,1)-x(1,3)
      xj12 = x(2,1)-x(2,3)
      xj21 = x(1,2)-x(1,3)
      xj22 = x(2,2)-x(2,3)
      det  = xj11*xj22-xj12*xj21
      raio = 1/3.d0*(x(1,1)+x(1,2)+x(1,3))
      if (det .le. 0.d0) goto 1000
c
c     Inversa da matriz Jacobiana:
c
      xji11 =  xj22/det
      xji12 = -xj12/det
      xji21 = -xj21/det
      xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao em relacao a r e z:
c       hx(i) = dN_i/dr,  hy(i) = dN_i/dz
c
      hx(1) =  xji11
      hx(2) =  xji12
      hx(3) = -xji11-xji12
      hy(1) =  xji21
      hy(2) =  xji22
      hy(3) = -xji21-xji22
c
c     Matriz constitutiva D (problema axissimetrico):
c
      my  = e(1)
      nu  = e(2)
      a   = 1.d0+nu
      b   = a*(1.d0-2.d0*nu)
      c   = my*(1.d0-nu)/b
      d11 = c
      d12 = my*nu/b
      d21 = d12
      d22 = c
      d33 = my/(2.d0*a)
c
c     Matriz de rigidez: Ke = 2 * pi * B^T D B * r * det/2
c       wt = Area do triangulo = det/2
c
      wt = 2.d0*PI*0.5d0*raio*det
      do 220 j = 1, 3
        k = (j-1)*2+1
        do 210 i = 1, 3
          l = (i-1)*2+1
c         s(l,k)     += (dNi/dx*d11*dNj/dx + dNi/dy*d33*dNj/dy)*A
          s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
c         s(l,k+1)   += (dNi/dx*d12*dNj/dy + dNi/dy*d33*dNj/dx)*A
          s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
c         s(l+1,k)   += (dNi/dy*d21*dNj/dx + dNi/dx*d33*dNj/dy)*A
          s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
c         s(l+1,k+1) += (dNi/dy*d22*dNj/dy + dNi/dx*d33*dNj/dx)*A
          s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
  210   continue
  220 continue
c
c     Vetor de forcas internas: p = Ke * u
c
      call lku(s,u,p,nst)
      return
 1000 continue
      print*, '*** ERRO elmt01: det nulo ou negativo no elemento', nel
      print*, '    Verifique a ordem dos nos (deve ser anti-horaria).'
      stop
      end

c ======================================================================
      subroutine actcol(a,b,jdiag,neq,afac,back)
c ======================================================================
c     Solver direto: decomposicao LtDL com armazenamento skyline.
c     Valido somente para matrizes simetricas.
c
c     Parametros de entrada:
c       a     = coeficientes da matriz (armazenados por altura de coluna)
c       b     = vetor independente
c       jdiag = vetor apontador do armazenamento skyline
c       neq   = numero de equacoes
c       afac  = .true. => fatoriza a matriz
c       back  = .true. => retrosubstitui
c
c     Parametros de saida:
c       a     = coeficientes da matriz triangularizada
c       b     = vetor solucao
c ======================================================================
      implicit real*8 (a-h,o-z)
      common/engys/ aengy
      dimension a(*),b(*),jdiag(*)
      logical afac,back
      aengy = 0.0d0
      jr = 0
      do 600 j = 1,neq
        jd = jdiag(j)
        jh = jd - jr
        is = j - jh + 2
        if(jh-2) 600,300,100
  100   if(.not.afac) goto 500
        ie = j - 1
        k  = jr + 2
        id = jdiag(is - 1)
        do 200 i = is, ie
          ir = id
          id = jdiag(i)
          ih = min0(id-ir-1,i-is+1)
          if(ih.gt.0) a(k) = a(k) - dot(a(k-ih),a(id-ih),ih)
  200   k = k + 1
  300   if(.not.afac) goto 500
        ir = jr+1
        ie = jd - 1
        k  = j - jd
        do 400 i = ir, ie
          id = jdiag(k+i)
          if(a(id).eq.0.0d0) goto 400
          d    = a(i)
          a(i) = a(i)/a(id)
          a(jd) = a(jd) - d*a(i)
  400   continue
  500   if(back) b(j) = b(j) - dot(a(jr+1),b(is-1),jh-1)
  600 jr = jd
      if(.not.back) return
      do 700 i = 1,neq
        id = jdiag(i)
        if(a(id).ne.0.0d0) b(i) = b(i)/a(id)
  700 aengy = aengy + b(i)*b(i)*a(id)
      j  = neq
      jd = jdiag(j)
  800 d  = b(j)
      j  = j - 1
      if(j.le.0) return
      jr = jdiag(j)
      if(jd-jr.le.1) goto 1000
      is = j - jd + jr + 2
      k  = jr - is + 1
      do 900 i = is,j
  900 b(i) = b(i) - a(i+k)*d
 1000 jd = jr
      goto 800
      end

c ======================================================================
      subroutine lku(s,u,p,nst)
c     Produto matriz-vetor local: p = s * u
c     Calcula o vetor de forcas internas do elemento.
c ======================================================================
      integer nst
      real*8 s(nst,nst),u(nst),p(nst)
      do 200 i = 1, nst
        p(i) = 0.d0
        do 100 j = 1, nst
          p(i) = p(i) + s(i,j) * u(j)
  100   continue
  200 continue
      return
      end

c ======================================================================
      real*8 function dot(a,b,n)
c     Produto interno (escalar) de dois vetores: dot = a . b
c ======================================================================
      integer n
      real*8 a(*), b(*)
      dot = 0.d0
      do 100 i = 1, n
        dot = dot + a(i) * b(i)
  100 continue
      return
      end

c ======================================================================
      subroutine mem(npos)
c     Verifica se ha memoria suficiente no vetor global.
c     Para se a posicao requerida ultrapassar o limite maximo.
c ======================================================================
      common /size/ max
      integer npos
      if ( (npos-1) .gt. max ) then
        print*, '*** ERRO: Memoria insuficiente!'
        print*, '    Posicao requerida:', npos-1,
     .          '  Maximo disponivel:', max
        print*, '    Aumente o parametro npos no programa principal.'
        stop
      endif
      return
      end

c ======================================================================
      subroutine mzero(m,i1,i2)
c     Zeragem de vetor inteiro: m(i1:i2) = 0
c ======================================================================
      integer i1,i2,m(*)
      do 100 i = 1, i2-i1+1
        m(i) = 0
  100 continue
      return
      end

c ======================================================================
      subroutine azero(a,i1,i2)
c     Zeragem de vetor real*8: a(i1:i2) = 0.d0
c ======================================================================
      integer i1,i2
      real*8 a(*)
      do 100 i = 1, i2-i1+1
        a(i) = 0.d0
  100 continue
      return
      end

c ======================================================================
      subroutine wdata(ix,id,x,f,u,nnode,numel,ndm,nen,ndf,nout)
c     Escrita dos resultados no arquivo de saida.
c     Formato compativel com visualizadores de pos-processamento.
c
c     Blocos gravados:
c       coor  => coordenadas nodais
c       elem  => conectividade dos elementos
c       nosc  => deslocamentos nodais
c       nvec  => campo vetorial com z=0 (para visualizacao 3D)
c       end   => marcador de fim de arquivo
c ======================================================================
      integer nnode,numel,ndm,nen,ndf,ix(nen+1,*),id(ndf,*)
      real*8  x(ndm,*),f(ndf,*),u(*),aux(6)
      write(nout,'(a,i5)') 'coor ', nnode
      do 100 i = 1, nnode
        write(nout,'(i10,3e15.5)') i,(x(j,i),j=1,2),0.0
  100 continue
      write(nout,'(a,i10)') 'elem ', numel
      do 200 i = 1, numel
        write(nout,'(10i10)') i,nen,(ix(j,i),j=1,nen)
  200 continue
c
c     Deslocamentos nodais:
c       Nos livres  (id>0): deslocamento calculado u(id)
c       Nos restritos (id=0): deslocamento prescrito f(j,i)
c
      write(nout,'(a,i2)') 'nosc ',ndf
      do 300 i = 1, nnode
        do 310 j = 1, ndf
          aux(j) = f(j,i)
          k = id(j,i)
          if(k .gt. 0) aux(j) = u(k)
  310   continue
        write(nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf)
  300 continue
c
c     Campo vetorial (adiciona z=0 para pos-processamento 3D):
c
      write(nout,'(a)') 'nvec '
      do 400 i = 1, nnode
        do 410 j = 1, ndf
          aux(j) = f(j,i)
          k = id(j,i)
          if(k .gt. 0) aux(j) = u(k)
  410   continue
        write(nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf),0.0d0
  400 continue
      write(nout,'(a)') 'end '
      return
      end

c ======================================================================
      subroutine pcg(a,b,jdiag,x,r,z,prec,d,neq)
c     Solver iterativo: Gradiente Conjugado Precondicionado (PCG)
c     Precondicionador: diagonal de K (precondicionador de Jacobi)
c
c     Resolve: K * x = b
c
c     Parametros:
c       a    = matriz de rigidez (formato skyline)
c       b    = vetor de forcas (entrada) / deslocamentos (saida)
c       jdiag= ponteiros do skyline
c       x    = vetor solucao (trabalho interno)
c       r    = vetor residuo (trabalho)
c       z    = vetor precondicionado (trabalho)
c       prec = precondicionador diagonal M = diag(K)
c       d    = direcao de descida conjugada (trabalho)
c       neq  = numero de equacoes
c ======================================================================
      implicit none
      integer neq,maxit,i,iter
      integer jdiag(*)
      real*8  a(*),prec(*),x(*),r(*),z(*),b(*),d(*)
      real*8  dnew,dold,alpha,beta,tol,norma_energia,dot
      external dot
      tol   = 1.0d-07
      maxit = 1000
c
c     Precondicionador diagonal: prec(i) = K(i,i) = a(jdiag(i))
c
      do 10 i = 1, neq
        prec(i) = a(jdiag(i))
        if (dabs(prec(i)) .lt. 1.d-30) then
          print*, '*** AVISO pcg: diagonal nula na equacao', i
          print*, '    Verifique as condicoes de contorno.'
          prec(i) = 1.d-30
        endif
   10 continue
c
c     Inicializacao: x=0, r=b-K*0=b, d=M^-1*r
c
      do 20 i = 1, neq
        x(i) = 0.d0
   20 continue
      call matvec(neq,a,jdiag,x,z)
      do 30 i = 1, neq
        r(i) = b(i) - z(i)
        d(i) = r(i) / prec(i)
   30 continue
      dnew = dot(r,d,neq)
      tol  = tol*tol*(dabs(dnew))
c
c     Loop principal do PCG:
c
      do 200 iter = 1, maxit
        call matvec(neq,a,jdiag,d,z)
        alpha = dnew / dot(d,z,neq)
        do 100 i = 1, neq
          x(i) = x(i) + alpha * d(i)
          r(i) = r(i) - alpha * z(i)
          z(i) = r(i) / prec(i)
  100   continue
        dold = dnew
        beta = dot(r,z,neq) / dold
        dnew = dot(r,z,neq)
        do 150 i = 1, neq
          d(i) = z(i) + beta * d(i)
  150   continue
c
c       Criterio de convergencia: ||r||_M < tol
c
        if (dabs(dnew) .lt. tol) then
          print*, 'PCG convergiu em', iter, 'iteracoes.'
          goto 300
        endif
  200 continue
c
c     Nao convergiu:
c
      print*, '*** AVISO pcg: nao convergiu em', maxit, 'iteracoes.'
      print*, '    Residuo final (dnew):', dnew
      print*, '    Tolerancia  (tol)   :', tol
c
c     Norma de energia: ||u||_K = sqrt(u^T K u)
c
  300 continue
      call matvec(neq,a,jdiag,x,z)
      norma_energia = dsqrt(dabs(dot(x,z,neq)))
      write(*,'(a,i6)')    'Numero de equacoes    :', neq
      write(*,'(a,f15.8)') 'Norma de energia      :', norma_energia
c
c     Copia solucao de volta para b (vetor de deslocamentos):
c
      do 400 i = 1, neq
        b(i) = x(i)
  400 continue
      return
      end

c ======================================================================
      subroutine matvec(neq,a,jdiag,x,y)
c     Produto matriz-vetor para matriz simetrica no formato skyline.
c     Calcula: y = K * x
c
c     Aproveita a simetria: percorre apenas a parte armazenada
c     (triangulo inferior de cada coluna) e contribui tanto para
c     y(i) quanto para y(i-l), explorando K(i,i-l) = K(i-l,i).
c ======================================================================
      implicit none
      integer neq,i,j,l,m
      integer jdiag(*)
      real*8  a(*),x(*),y(*)
      do 50 i = 1, neq
        y(i) = 0.d0
   50 continue
      do 200 i = 1, neq
c
c       Produto diagonal: y(i) += K(i,i) * x(i)
c
        j    = jdiag(i)
        y(i) = a(j)*x(i)
c
c       Produtos fora da diagonal (explorando simetria de K):
c       K(i, i-l) = a(j-l) contribui para y(i) e y(i-l)
c
        if (i .le. 1) goto 200
        m = jdiag(i)-jdiag(i-1)-1
        do 100 l = 1, m
          y(i-l) = y(i-l) + a(j-l)*x(i)
          y(i)   = y(i)   + a(j-l)*x(i-l)
  100   continue
  200 continue
      return
      end

c ======================================================================
      subroutine wvtk(ix,id,x,f,u,sxa,sya,txya,nnode,numel,ndm,nen,
     .                ndf,nvtk)
c     Escrita dos resultados em formato VTK para visualizacao no ParaView.
c
c     Formato: VTK Unstructured Grid ASCII
c       - Coordenadas nodais (com z=0 para 2D)
c       - Conectividade dos elementos
c       - Dados nodais: deslocamentos (u_x, u_y, u_z)
c       - Dados por elemento: tensoes (sigma_x, sigma_y, tau_xy)
c ======================================================================
      integer nnode,numel,ndm,nen,ndf,ix(nen+1,*),id(ndf,*),nvtk
      real*8  x(ndm,*),f(ndf,*),u(*),aux(3)
      real*8  sxa(*),sya(*),txya(*)
      integer i,j,k
c
c     Cabecalho VTK:
c
      write(nvtk,'(a)') '# vtk DataFile Version 2.0'
      write(nvtk,'(a)') 'MEF - Elasticidade Plana 2D'
      write(nvtk,'(a)') 'ASCII'
      write(nvtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
c
c     Bloco de coordenadas nodais:
c
      write(nvtk,'(a,i8,a)') 'POINTS ', nnode, ' float'
      do 100 i = 1, nnode
        write(nvtk,'(3e16.8)') x(1,i), x(2,i), 0.0d0
  100 continue
c
c     Bloco de conectividade dos elementos:
c     Total de valores = numel*(nen+1)  [numero de nos + tipo]
c     Tipo VTK: 5 para triangulos (T3), 9 para quadrilateros (Q4)
c
      write(nvtk,'(a,2i8)') 'CELLS ', numel, numel*(nen+1)
      do 200 i = 1, numel
        write(nvtk,'(10i8)') nen, (ix(j,i)-1, j=1,nen)
  200 continue
c
c     Tipos de celulas VTK:
c     5 = VTK_TRIANGLE (T3)
c     9 = VTK_QUAD (Q4)
c
      write(nvtk,'(a,i8)') 'CELL_TYPES ', numel
      if (nen .eq. 3) then
        do 300 i = 1, numel
          write(nvtk,'(i2)') 5
  300   continue
      else if (nen .eq. 4) then
        do 350 i = 1, numel
          write(nvtk,'(i2)') 9
  350   continue
      endif
c
c     Dados por elemento: tensoes (sigma_x, sigma_y, tau_xy)
c
      if (ndf .eq. 2) then
        write(nvtk,'(a,i8)') 'CELL_DATA ', numel
        write(nvtk,'(a)') 'SCALARS sigma_x float 1'
        write(nvtk,'(a)') 'LOOKUP_TABLE default'
        do 600 i = 1, numel
          write(nvtk,'(1e16.8)') sxa(i)
  600   continue
        write(nvtk,'(a)') 'SCALARS sigma_y float 1'
        write(nvtk,'(a)') 'LOOKUP_TABLE default'
        do 610 i = 1, numel
          write(nvtk,'(1e16.8)') sya(i)
  610   continue
        write(nvtk,'(a)') 'SCALARS tau_xy float 1'
        write(nvtk,'(a)') 'LOOKUP_TABLE default'
        do 620 i = 1, numel
          write(nvtk,'(1e16.8)') txya(i)
  620   continue
      endif
      if (ndf .eq. 2) then
c
c     Dados nodais: deslocamentos (u_x, u_y, u_z com z=0)
c
         write(nvtk,'(a,i8)') 'POINT_DATA ', nnode
         write(nvtk,'(a)') 'VECTORS Deslocamentos float'
         do 400 i = 1, nnode
           aux(1) = f(1,i)
           aux(2) = f(2,i)
           aux(3) = 0.0d0
           k = id(1,i)
           if (k .gt. 0) aux(1) = u(k)
           k = id(2,i)
           if (k .gt. 0) aux(2) = u(k)
           write(nvtk,'(3e16.8)') aux(1), aux(2), aux(3)
  400    continue
         return
      elseif (ndf .eq. 1) then
c
c     Dados nodais: temperatura
c
         write(nvtk,'(a,i8)') 'POINT_DATA ', nnode
         write(nvtk,'(a)') 'SCALARS Temperatura float'
         write(nvtk,'(a)') 'LOOKUP_TABLE default'
         do 500 i = 1, nnode
           aux(1) = f(1,i)
           k = id(1,i)
           if (k .gt. 0) aux(1) = u(k)
           write(nvtk,'(1e16.8)') aux(1)
  500    continue
         return
      endif
      end