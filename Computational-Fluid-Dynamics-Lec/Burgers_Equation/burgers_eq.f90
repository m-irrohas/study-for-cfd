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
    integer num, n_step, i, j
    real x_min, x_max, u1, u2, cfl, dx, dt , x_i, c, lx, dtdx
    ! 格子の最大点と最小点
    x_min = -5.0
    x_max = 5.0
    num = 100 !刻み数

    ! その他
    n_step = 50
    dx = (x_max-x_min)/num
    print '(I5)', n_step
    
    ! 格子
    x_i = x_min
    x(1) = x_i
    do i=1, num
        x_i = x_i + dx
        x(i) = x_i
    end do

    u1 = 1.0
    u2 = -1.0
    call solve_godunov(u1, u2, "./output/godunov1.csv")
    call solve_engquist_osher(u1, u2, "./output/engquist_osher1.csv")

    u1 = -1.0
    u2 = 1.0
    call solve_godunov(u1, u2, "./output/godunov2.csv")
    call solve_engquist_osher(u1, u2, "./output/engquist_osher2.csv")

contains

    real function godunov_func(u, u_next)
        real u, u_next, c, a
        c = (u+u_next)/2.
        if (u<0 .and. u_next<0) then
            a = 1.0/2.0 *u_next**2
        else if (u>=0 .and. u_next>=0) then
            a = 1.0/2.0 *u**2
        else if (u>=0 .and. 0>u_next .and. c<0) then
            a = 1.0/2.0 *u_next**2
        else if (u>=0 .and. 0>u_next .and. c>=0) then
            a = 1.0/2.0 *u**2
        else
            a = 0
        end if
        godunov_func = a
    end function godunov_func

    real function engquist_osher_func(u, u_next)
        real u, u_next, c, a
        c = (u+u_next)/2.
        if (u<0 .and. u_next<0) then
            a = 1.0/2.0 *u_next**2
        else if (u>=0 .and. u_next>=0) then
            a = 1.0/2.0 *u**2
        else if (u>=0 .and. 0>u_next) then
            a = 1.0/2.0 *(u**2+u_next**2)
        else
            a = 0
        end if
    engquist_osher_func = a
    end function engquist_osher_func

    subroutine solve_godunov(u1_init, u2_init, output)
        ! ゴドノフ法を解く
        ! Arg:
        !   u1_init(float) 初期条件(x<=0)
        !   u2_init(float) 初期条件(x>0)
        ! Return:
        !   csv => output
        !       row[0]は位置
        !       row[1:]は速度の遷移
        real u1_init
        real u2_init
        real f_plus, f_minus
        real,dimension(100)::u_before
        character(*) output
        ! 初期条件
        do i=1, num
            if(x(i)<=0) then
                u_init(i) = u1_init
            else
                u_init(i) = u2_init
            end if
            if (-1E-4<x(i) .and. x(i)<1E-4) then
                u_init(i) = 0
            end if
        end do
        !初期条件をコピー
        do i=1, num
            u(i) = u_init(i)
        end do

        !出力ファイル
        open(20, file=output, status='replace')
        write (20,*) (x(i), i=1,num)
        write(20,*) (u_init(i), i=1,num)
        do j=1, n_step !時間ステップ
            do i = 1, num
                u_before(i) = u(i) !演算で使うのは今の状態
            end do
            do i=2,num-1 !位置
                f_plus = godunov_func(u_before(i), u_before(i+1))
                f_minus = godunov_func(u_before(i-1), u_before(i))
                u(i) = u_before(i) - (f_plus-f_minus)
            end do
            write (20,*) (u(i),i=1,num)
        end do
        close(20)
    end subroutine

    subroutine solve_engquist_osher(u1_init, u2_init, output)
        ! Engquist-Osher法を解く
        ! Arg:
        !   u1_init(float) 初期条件(x<=0)
        !   u2_init(float) 初期条件(x>0)
        ! Return:
        !   csv => output
        !       row[0]は位置
        !       row[1:]は速度の遷移
        real u1_init
        real u2_init
        real f_plus, f_minus
        real,dimension(100)::u_before
        character(*) output
        ! 初期条件
        do i=1, num
            if(x(i)<=0) then
                u_init(i) = u1_init
            else
                u_init(i) = u2_init
            end if
        end do
        !初期条件をコピー
        do i=1, num
            u(i) = u_init(i)
        end do

        !出力ファイル
        open(18, file=output, status='replace')
        write (18,*) (x(i), i=1,num)
        write(18,*) (u_init(i), i=1,num)
        do j=1, n_step !時間ステップ
            do i = 1, num
                u_before(i) = u(i) !演算で使うのは今の状態
            end do
            do i=2,num-1 !位置
                f_plus = engquist_osher_func(u_before(i), u_before(i+1))
                f_minus = engquist_osher_func(u_before(i-1), u_before(i))
                u(i) = u_before(i) - (f_plus-f_minus)
            end do
            write (18,*) (u(i),i=1,num)
        end do
        close(18)
    end subroutine
end program main