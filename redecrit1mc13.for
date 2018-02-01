c     programa redecrit1mc13 
c     igual a redecrit1mc11, mas usa a versao madchar13 para o calculo da matriz de vizinhanca 
c     calcula distancia entre duas redes determinadas para valores distintos 
c     de uma variável de controle 
c     a partir de suas matrizes de vizinhanca calculada com madchar13 
c     entra dados no formato de uma matriz de interação dependente do periodo 
c     unica diferença é escrever a matriz de vizinhança para o último valor de ip
c     entra com variaveis reais, ao passo que redecrit0mc11 entra variaveis inteiras

      parameter(npm=10000)
      integer*1 am(npm,npm)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2),kki(0:npm,0:npm)
      real lar(2,npm)
      real idis(npm,npm)
      character*100 entrada3,saida5,saida6,saida7

      open (unit=2,file='redecrit1mc13.dat')

c     nm:numero maximo de nos na rede 
c     np:potencia booleana maxima 
c     pmin,delp,ns: valor minimo, intervalo de variacao do parametro, numero de valores do parametro
c     nc:numero de nos conectados na rede 
c     na:numero de arestas na rede

c     entrada de dados
c==================================================================

      read (2,*)nm,np
      read (2,*)pmin,delp,ns
      read (2,500)entrada3
      read (2,500)saida5
      read (2,500)saida6
      read (2,500)saida7

      open (unit=3,file=entrada3)
      open (unit=5,file=saida5)
      open (unit=6,file=saida6)
      open (unit=7,file=saida7)

      write(5,*)'      rede   nc   na  mcl nmcl diam        xmed    d(i,
     .i+1)'
      write(6,*)'      rede   nc   na  mcl nmcl diam        xmed    d(i,
     .i+1)'

      do i = 1,npm
       lar(1,i) = 0
      enddo

      do i = 1,nm
       read (3,*)(lar(1,j),j=1,nm)
       do j = 1,nm
        idis(i,j) = lar(1,j)
       enddo
      enddo

      xp = pmin
      na = 0
      nc = 0
      do i = 1,nm
       ic = 0
       do j = 1,nm
        am(i,j) = 0
c        if(idis(i,j).le.xp.and.idis(i,j).gt.0)then 
        if(idis(i,j).ge.xp.and.idis(i,j).gt.0)then 
         am(i,j) = 1
         na = na + 1
         ic = 1
        endif
       enddo
       nc = nc + ic
      enddo
      na = na/2

      call madchar13(am,mv1,nm,nc,np,kki,xmed,id,mcl,idege)

      is = 1
      d = 0.
      write(5,520)xp,nc,na,mcl,idege,id,xmed,d/delp
      write(*,520)xp,nc,na,mcl,idege,id,xmed,d/delp

c================================================================== 
c     comeco do grand loop no numero de amostras
c==================================================================

      do is = 1,ns

       xp = pmin + is*delp

       na = 0
       nc = 0
       do i = 1,nm
        ic = 0
        do j = 1,nm
         am(i,j) = 0
c        if(idis(i,j).le.xp.and.idis(i,j).gt.0)then 
        if(idis(i,j).ge.xp.and.idis(i,j).gt.0)then 
          am(i,j) = 1
          na = na + 1
          ic = 1
         endif
        enddo
        nc = nc + ic
       enddo
       na = na/2

       call madchar13(am,mv2,nm,nc,np,kki,xmed,id,mcl,idege)

       d = dist(mv1,mv2,nm)

       write(5,520)xp,nc,na,mcl,idege,id,xmed,d/delp
       write(*,520)xp,nc,na,mcl,idege,id,xmed,d/delp

       do i = 1,nm*(nm-1)/2
         mv1(i) = mv2(i)
       enddo
      enddo

      ll = 0
      do i = 1,nm
       do j = i+1,nm
        ll = ll + 1
        idis(i,j) = mv1(ll)
        idis(j,i) = mv1(ll)
       enddo
      enddo

      do i = 1,nm
       write(6,540)(int(idis(i,j)),j=1,nm)
      enddo

      do i = 1,nm
       do j = 1,nm
        if(idis(i,j).ne.1)am(i,j)=0
       enddo
      enddo

      do i = 1,nm
       write(7,510)(am(i,j),j=1,nm)
      enddo

500   format(a100) 
510   format(10000i1) 
520	format(1x,f10.4,5(1x,i7),2(2x,e10.3)) 
c520	format(1x,f10.8,5(1x,i4),2(2x,e10.3)) 
c520  format(6(1x,i4),2(2x,e10.3))
c530   format(10000(f5.0,1x)) 
c530   format(10000(i6,1x)) 
540	format(10000i4)

      stop

      end

c     fim do programa principal
c==================================================================
c==================================================================

      subroutine madchar13(a,mv,nm,nc,np,kk,xlmd,id,mcl,idege)
c==================================================================

      parameter(npm=10000)
      integer*1 a(npm,npm)
      integer*2 mv(npm*(npm-1)/2),kk(0:npm,0:npm),lla(npm),pro(npm)
      integer*2 pra(npm)
      real lis(80*npm,2)
      real ga(0:npm), kmd(0:npm)
      real gamd(0:npm),lmd(0:npm),set(npm,npm)
c================================================================== 
c     coloca zero no valores medios das distancias entre nos
      xlmd = 0.
      ylmd = 0.
      zlmd = 0.

      do i = 1,npm
       lla(i) = 0
       do j = 1,npm
        set(i,j) = 0.
       enddo
      enddo

      do i = 1,(npm*npm-1)/2
       mv(i) = 0
      enddo

      do i = 0,npm
       lmd(i) = 0.
       ga(i) = 0
       kmd(i) = 0
       gamd(i) = 0
       do j = 0,npm
        kk(i,j) = 0
       enddo
      enddo
c================================================================== 
c     escreve a(i,j) em entrada

      do i = 1,nm-1
       do j = i+1,nm
        a(j,i) = a(i,j)
        ll = (2*nm-i)*(i-1)/2+j-i
        mv(ll) = a(i,j)
       enddo
      enddo
      do i = 1,nm
        a(i,i) = 1
        pra(i) = 1
        pro(i) = 0
      enddo
c================================================================== 
c     prepara a lista de elementos nao nulos de a(i,j)

      do i = 1,nm
       do j = 1,2
        lis(i,j) = 0.
       enddo
      enddo
      iq = 1
      do i = 1,nm
       lis(i,2) = iq
       do j = 1,nm
        if(a(i,j).eq.1)then
         lis(iq,1) = j
         iq = iq + 1
        endif
       enddo
      enddo
      lis(i,2) = iq
c================================================================== 
c     comeca determinacao da matriz vizinhanca em tmad1.dat
c================================================================== 
c     soma matriz identidade com mad e coloca matriz com todos os vizinhos em tmad5
c================================================================== 
c     comeca o loop para calculo as propriedades das diferentes matrizes mad(ip)

      do i = 1,nm-1
       do ip = 1,np
        if(pra(i).eq.0.and.pro(i).lt.ip)goto 90 
c================================================================== 
c	calcula coeficiente de clusterização e grau da matriz mad(ip)

C        call rede1(a,nm,npm,kk,ga,ip)
C
C        if (ip.eq.1) then
C         gamd(ip) = ga(0)/(0*nm+nc)
C         do i = 1,nm
C          set(i,ip) = ga(i)
C         enddo
C        endif
C        kmd(ip) = float(kk(0,ip))/(0*nm+nc)
C        do i = 1,nm
C         lmd(i)= lmd(i) + kk(i,ip)*ip
C        enddo
C        kk(0,ip) = 0
C        do i = 1,nm
C         kk(0,ip) = kk(0,ip) + pra(i)
C        enddo
c================================================================== 
c	calcula grau de assortatividadeda rk matriz mad(ip)

c        call assorta(am,nm,npm,nvm,kk,ip,rrk) 
c        rk(ip) = rrk
c================================================================== 
c	desvio para finalizar o programa

C        if (kk(0,ip).eq.0) go to 100
        if (ip.eq.np) go to 100
c================================================================== 
c	comeca a obtencao de am**ip - le mad e transfere para am
c================================================================== 
c	le a matriz tmad5 faz pb tmad5*mad e guarda resultado em tmad6
c==================================================================

        pra(i) = 0

        do j = i+1,nm

         a(i,j) = 0

         if (a(j,i).eq.1) goto 126

         do k = lis(j,2),lis(j+1,2)-1
	    k1 = lis(k,1)
          if(i.lt.k1) then 
           if (mv((2*nm-i)*(i-1)/2+k1-i).gt.0)then
	      if(mv((2*nm-i)*(i-1)/2+k1-i).le.ip)go to 124
           endif
          else if (i.gt.k1) then 
           if (mv((2*nm-k1)*(k1-1)/2+i-k1).gt.0)then
            if(mv((2*nm-k1)*(k1-1)/2+i-k1).le.ip)go to 124
           endif
          endif
         enddo

         goto 125
124      a(i,j) = 1
         pra(i) = 1
125      continue
         ll = (2*nm-i)*(i-1)/2+j-i
         mv(ll) = a(i,j)*(ip+1)
         pro(j) = max(pro(j),mv(ll))

126     enddo

c================================================================== 
c     obtem a nova matriz com vizinhos exclusivos de ordem ip+1: am=tmad6-tmad5
c================================================================== 
c     transfere a nova matriz com todos os vizinhos de tdma1 para tdma0
c================================================================== 
c     atualiza a matriz de adjacencia ponderada, le tmad1, soma (ip+1)*am e coloca em tmad2
c================================================================== 
c     transfere a nova matriz com todos os vizinhos ponderados de tmad2 para tmad1 
c==================================================================

        do j = i+1,nm
         a(j,i) = a(j,i) + a(i,j)
        enddo

       enddo

c====================================================== 
c     fim do loop em ip

90    enddo
c====================================================== 
c     saida dos resultados para a amostra is

100   npp = 1

      do i = 1,nm-1
       do j = i+1,nm
        ll = (2*nm-i)*(i-1)/2+j-i
        npp = max(npp,mv(ll))
       enddo
      enddo

c======================================================
c     calcula a ordem de sub-grafo de cada no

      do i = 1,nm
c       a(i,i) = 0
       do j = i+1,nm
        a(i,j) = a(j,i)
       enddo
      enddo


      do i = 1,nm
       kk(i,0) = 0
       do j = 1,nm
        kk(i,0) = kk(i,0) + a(i,j)
       enddo
      enddo

      do i = 1,nm
c       a(i,i) = 0
       kk(i,i) = 0
      enddo

      ll = 0
      do i = 1,nm-1
       do j = i+1,nm
        ll = ll + 1
c        a(i,j) = mv(ll)
c        a(j,i) = mv(ll)
        kk(i,j) = mv(ll)
        kk(j,i) = mv(ll)
       enddo
      enddo

      do i = 1,nm
       lmd(i) = 0
       do j = 1,nm
c        lmd(i) = lmd(i) + a(i,j)
        lmd(i) = lmd(i) + kk(i,j)
       enddo	
      enddo	
c====================================================== 
c     calcula a numero de nos no maior cluster

      mcl = 0
      do i = 1,nm
       mcl = max(kk(i,0),mcl)
      enddo
c====================================================== 
c     calcula a distancia minima média da rede limitada ao(s) maior(es) cluster(s)

c     nm2 = (0*nm+nc)*(0*nm+nc)
      xm2 = 0.
      if (mcl.gt.1)xm2 = 1./float((mcl-1)*mcl)
      imc = 0
      do i = 1,nm
       if (kk(i,0).eq.mcl)then
        imc = imc + 1
        xlmd = xlmd + lmd(i) 
c       xlmd = xlmd + lmd(i)*xm2
        ylmd = ylmd + lmd(i)/kk(i,0)/kk(i,0)
        zlmd = zlmd + lmd(i)/kk(i,0)/(0*nm+nc)
       endif
      enddo
      idege = imc/mcl
c      write(5,*)xlmd,xm2,mcl,imc
      xlmd = xlmd*xm2*mcl/imc
c      write(5,*)xlmd,xm2,mcl,imc
c====================================================== 
c     transfere a matriz ponderada para o campo mv(np,np)

c====================================================== 
c     calcula diametro

      id = 0
      do i = 1,nm*(nm-1)/2
       id = max(mv(i),id)
      enddo
c====================================================== 
c     calcula dimensao fractal 
c====================================================== 
c     escreve distancia minima media de cada no em saida9
c====================================================== 
c     escreve coeficiente de clusterização em saida7
c====================================================== 
c     escreve coeficientes de assortatividade em saida8
c====================================================== 
c     escreve distribuição de nos em saida8
c=====================================================
      return

      end

c     fim da subroutine madchar13(a,mv,nm,np,xmd,id)
c=====================================================================
c=====================================================================

      subroutine rede1(a,n,np,kk,ga,ip)
c=====================================================================
      parameter(npm=10000)
      integer*1 a(npm,npm)
      integer*2 kk(0:npm,0:npm)
      real ga(0:npm)
      integer vz(np)
c===================================================================== 
c     calcula o grau de cada no

      kk(0,ip) = 0
      do i = 1,n
       kk(i,ip) = 0
       do j = 1,n
        if(i.ne.j)then
         ii = min(i,j)
         jj = max(i,j)
         kk(i,ip) = kk(i,ip) + a(ii,jj)
        endif
       enddo
       kk(0,ip) = kk(0,ip) + kk(i,ip)
      enddo

      if (ip.gt.0) return

c===================================================================== 
c     calcula o coeficiente de aglomeração de cada noh

      if (ip.gt.1)return

      ga(0) = 0.
      do i=1,n
       do j=1,n
        vz(j) = -1
       enddo

       ordem = 0
       do j=1,n
        if (a(i,j) .eq. 1) then
         ordem=ordem+1
         vz(ordem) = j
        endif
       enddo

       ctotal = ((ordem * (ordem-1))/2)
       sumg=0
       do j=1,ordem
        li = vz(j)
        do k=1,n
         if (a(li,k).eq.1.and.a(i,k).eq.1.and.li.ne.k.and.k.ne.i)then
          sumg=sumg+1
         endif
        enddo
       enddo

       clocal = (sumg)/2
       if (ctotal .gt. 0) then
        ga(i) = clocal/ctotal
       else
        ga(i)=0
       endif
       ga(0) = ga(0) + ga(i)
      enddo

      return
      end

c     fim da subroutine rede1(a,n,np,nvm,kk,ga,ip)
c=====================================================================
      real function dist(mv1,mv2,nm)

      parameter(npm=10000)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2)

      dist = 0.

      do i = 1,nm*(nm-1)/2
       dist = dist + 2*(mv1(i)-mv2(i))**2
      enddo

      dist = sqrt(dist)/(nm-1)/nm

      return
      end

c     fim da real function dist(mv1,mv2,nm)
c=====================================================================
c=====================================================================
