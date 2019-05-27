program main
    implicit none
    ! 格子定義
    real,dimension(100)::x
    real,dimension(100)::u
    real,dimension(100)::u_init
    real,dimension(100)::u_
    integer num, nlast, i, j
    real x_min, x_max, u1, u2, a, cfl, dx, dt ,l_x, x_i
    ! 格子の最大点と最小点
    x_min = -5.0
    x_max = 5.0
    num = 100

    ! 速度の定義
    u1 = 1.0
    u2 = 0.0

    ! その他
    a = 1.0
    nlast = 10
    cfl = 0.5
    dx = (x_max-x_min)/num
    dt = 0.2
    l_x = dt*a*nlast
    ! 格子
    x_i = x_min
    x(1) = x_i
    do i=1, num
        x_i = x_i + dx
        x(i) = x_i
    end do

    ! 初期条件
    do i=1, num
        if(x(i)<=0) then
            u_init(i) = u1
        else
            u_init(i) = u2
        end if
    end do

    !!!FTCS
    !初期条件をコピー
    do i=1, num
        u(i) = u_init(i)
    end do

    !出力ファイル
    open(18, file='test.csv', status='replace')
    write (18,*) (x(i), i=1,num)
    do j=1, nlast !時間ステップ
        do i=2,num-1 !位置
            u(i) = u(i)-cfl/2*(u(i+1)-u(i-1))
        end do
        write (18,*) (u(i),i=1,num)
    end do
    close(18)


end program main
! 格子生成
