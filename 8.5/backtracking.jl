function strong_backtracking(f, ∇f, x0, d; α0=1.0, β=1e-4, σ=0.9)
    y0, g0 = f(x0), ∇f(x0)⋅d
    @assert g0 < 0 "Направление d не убывающее"

    α = α0
    yprev, αprev = y0, 0.0
    αlo = αhi = NaN

    # bracketing phase
    # поиск интервала (`αlo`, `αhi`), в котором существует α,
    # удовлетворяющее строгим условиям Вольфе
    for i in 1:100
        y = f(x0 + α*d)
        # первое условие Вольфе нарушено или функция увеличилась
        if y > y0 + β*α*g0 || y ≥ yprev
            αlo, αhi = αprev, α
            break
        end
        g = ∇f(x0 + α*d)⋅d
        if abs(g) ≤ -σ*g0  # 2ое условие Вольфе выполнено, `α` найдено
            return α
        elseif g ≥ 0
            αlo, αhi = αprev, α
            break
        end
        yprev, αprev, α = y, α, 2α
    end
    # zoom phase
    # поиск α, для которого выполняются условия Вольфе
    ylo = f(x0 + αlo*d)
    for i in 1:100
        # бинарный поиск, но может быть и интерполяция
        α = (αlo + αhi)/2

        y = f(x0 + α*d)
        if y > y0 + β*α*g0 || y ≥ ylo  # `αhi` слишком большое
            αhi = α
        else
            g = ∇f(x0 + α*d)⋅d
            # 2ое условие Вольфе удовлетворено, `α` найдено
            abs(g) ≤ -σ*g0 && return α
            # обновление границ по бинарному поиску на производную
            if g ≥ 0
                αhi = α
            else
                αlo = α
            end
        end
    end
    @error "Ошибка в линейном поиске, α не найдено"
    return NaN
end