program main
    ! 2つの初期条件に対して，安定・不安定の計4通りの結果を出力する。
    ! Arg:
    !   None(内部で勝手にやる)
    ! Return:
    !   *.csv => output直下に
    !       shape=(step数, xの個数)
    implicit none
    ! 格子定義
    real,dimension(100)::x
    real,dimension(100)::u
    real,dimension(100)::u_init
    real,dimension(100)::u_
    integer num, n_step, i, j
    real x_min, x_max, u1, u2, a, cfl, dx, dt ,l_x, x_i, c, lx
    ! 格子の最大点と最小点
    x_min = -5.0
    x_max = 5.0
    num = 100 !刻み数

    ! その他
    cfl = 0.5 ! = c*(dt/dx)
    dx = (x_max-x_min)/num !刻みから逆算
    dt = 1
    c = cfl*dx/dt !クーラン数から逆算
    lx = 2.0
    n_step = int(lx/dt/c) !移動距離から逆算
    print '(I5)', n_step
    ! 格子
    x_i = x_min
    x(1) = x_i
    do i=1, num
        x_i = x_i + dx
        x(i) = x_i
    end do

    ! 速度の定義
    u1 = 1.0
    u2 = 0.0
    call solve_ftcs(u1, u2, "../../output/ftcs1.csv")

    u1 = 0.0
    u2 = 1.0
    call solve_ftcs(u1, u2, "../../output/ftcs2.csv")
contains

    subroutine solve_ftcs(u1_init, u2_init, output)
        ! FTCSを解く
        ! Arg:
        !   u1_init(float) 初期条件(x<=0)
        !   u2_init(float) 初期条件(x>0)
        ! Return:
        !   csv => output
        !       row[0]は位置
        !       row[1:]は速度の遷移
        real u1_init
        real u2_init
        character(*) output
        ! 初期条件
        do i=1, num
            if(x(i)<=0) then
                u_init(i) = u1_init
            else
                u_init(i) = u2_init
            end if
        end do

        !!!FTCS
        !初期条件をコピー
        do i=1, num
            u(i) = u_init(i)
        end do

        !出力ファイル
        open(20, file=output, status='replace')
        write (20,*) (x(i), i=1,num)
        write(20,*) (u_init(i), i=1,num)
        do j=1, n_step !時間ステップ
            do i=2,num-1 !位置
                u(i) = u(i)-cfl/2*(u(i+1)-u(i-1))
            end do
            write (20,*) (u(i),i=1,num)
        end do
        close(20)
    end subroutine

end program main