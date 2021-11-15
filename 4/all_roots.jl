using Polynomials

function ridders(f, x₁, x₂; maxiter=25, xtol=eps()*10, ftol=eps()*10)
	if x₁ > x₂; x₁, x₂ = x₂, x₁; end
	y₁, y₂ = f(x₁), f(x₂)
	y₁ * y₂ > 0 && return error("Функция должна иметь разные знаки в концах отрезка")
	y₁ == 0 && return x₁
	y₂ == 0 && return x₂
	
	for i in 1:maxiter
		xmid = (x₁ + x₂) / 2
		ymid = f(xmid)
		xnew = xmid + (xmid - x₁) * sign(y₁) * ymid / sqrt(ymid^2 - y₁*y₂)
		# xnew = xmid + (xmid - x₁) * (ymid / y₁) / sqrt((ymid / y₁)^2 - y₂/y₁)  # без функции sign
		ynew = f(xnew)
		
		ynew == 0 && return xnew
		
		if ynew * y₁ < 0
			x₂, y₂ = xnew, ynew
		elseif ynew * y₂ < 0
			x₁, y₁ = xnew, ynew
		end
		if abs(ynew) < ftol || abs(x₁ - x₂) < xtol
			return xnew
		end
	end
	return error("Число итераций превышено.")
end

# аргумент roots показывает, сколько корней нужно найти
function all_roots(f, a, b, roots)
	locations = zeros(roots)
	h = 0.1

	# локализация корней
	while true
		i = 1
		for x in range(a, b, step = h)
			if f(x)*f(x+h) <= 0
				if i > 1 && locations[i-1] == x-h
					continue
				end

				locations[i] = x
				i += 1
			end
			if i >= roots+1
				break
			end
		end

		if i >= roots+1
			break
		else
			h /= 2
		end

		if h < 0.0001
			error("не работает")
		end
	end

	# уточнение корней
	for i in 1:roots
		locations[i] = ridders(f, locations[i], locations[i]+h)
	end

	return locations
end