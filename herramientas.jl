__precompile__() # Se precompila el paquete.

module herramientas

export metodo_newton

"""documentación del método de Newton"""

function metodo_newton(F,df,valorinicial) #Se define una función que calcule las raíces 
    #usando el método de Newton que tome como entradas a una función cualquiera, la 
    #derivada y el valor inicial.
    x=valorinicial
    for i in 1:100
    x=x-F(x)/df(x);
    end
   return x #Se regresa el valor de la raiz.
end

export metodo_newton2

"""documentación del segundo método de Newton"""

function metodo_newton2(F,xi,h,error=1e-8) #Se define otra función que calcule las raices 
    #usando el método de Newton nuevamente, pero ahora sólo tomando las raices hasta 
    #cierto valor denotado por "error".
    df=F(xi+h)-F(xi)
    d=error+1
    while error<d
        xi=xi-F(xi)/(df)
        d=(F(xi)^2)^(1/2)
    end
    return xi #Regresa la raiz.
end

export metodo_newton3

"""documentación del tercer método de Newton"""

function metodo_newton3(F,xa,xb,h,error=1e-10) #Se crea una rutina donde en lugar de 
    #tomar una sóla condición inicial, reciba un rango de condiciones iniciales.
    raices=[] #Se define el vector raices.
    for i in linspace(xa,xb,200) #Se realiza la iteración 200 veces en el intervalo inicial.
        
        d=error+1
        while error<d
            i=i-F(i)/(F(i+h)-F(i))
            d=(F(i)^2)^(1/2)
        end
     
        push!(raices,i) #Guarda las raices en el vector definido raices.
    end
    return raices #Devuelve el vector.
end

export metodo_trapecio

"""documentación del método del Trapecio"""

function metodo_trapecio(f,a,b,n) #Se define una función que calcule la integral de una 
    #función dada usando el método del trapecio que tome como entradas a la función, el 
    #intervalo de integración, y el número de particiones del intervalo.
    x=a
    I=0
    while x<b
        x=x+((b-a)/n)
        I += (b-a)*(f(x)+f(x+(b-a)/n))/2n
    end
   return I #Se regresa el valor de la integral.
end

export metodo_rectangulo

"""documentación del método del Rectángulo"""


function metodo_rectangulo(F,a,b,n) #Se define una función que calcule la integral de 
    #una función dada usando el método del rectángulo que tome como entradas a la función, 
    #el intervalo de integración, y el número de particiones del intervalo.
    x=a
    I=0
    while x<b
        x=x+((b-a)/n)
        I += ((b-a)/n)*F(x+((b-a)/2n))
    end
   return I #Se regresa el valor de la integral.
end

export metodo_simpson

"""documentación del método de Simpson"""

function metodo_simpson(F,a,b,n) #Se define una función que calcule la integral de una función 
    #dada usando el método de Simpson que tome como entradas a la función, el intervalo de 
    #integración, y el número de particiones del intervalo.
    x = a
    I = 0
    while x < b
        x = x + ((b-a) / n)
        I += (b - a) * (F(x) + 4 * F(x + (b-a) / 2n) + F(x + (b-a)/n)) / 6n
    end
   return I #Se regresa el valor de la integral.
end

export metodo_euler

"""documentación del método de Euler vectorial"""

function metodo_euler_vec(f,list,x0) #Se crea una función que implementa el método de Euler vectorial. 
    #Dicha función toma como entradas la función, una lista y el valor inicial x0.
     x = x0 #Se asigna a x el valor inicial x0.
     h = list[2]-list[1] #Se escribe el paso de h.
     listx = [] #Se crea un vector vacío.
     for i in 2:length(list)
        t = i*h #Se calculan los valores de t.
        x = x + f(x,t)*h #Se realiza la operación del método de Euler.
        push!(listx,x) #Se guarda en listx los valores de x que se obtienen del método de Euler.
     end
     return listx #Se regresa la lista de x.
end

export metodo_implicito_euler

"""documentación del método implícito de Euler"""

function metodo_implicito_euler(f,x0,t0,tf,h) #Se crea una función que permite, mediante el método
    # implícito de Euler, obtener la solución aproximada de una ecuación diferencial. Dicha función 
    #toma como entradas la función que aparece en tal ecuación, la condición inicial, el valor inicial 
    #de t, el valor final de t, y el valor de h.
    listt = [] #Se crea un vector cuyo valor inicial es t0 y valor final tf donde la separación entre  
    #sus valores está dada por el valor de h.
    listx = [] #Se crea un vector en donde se guardan los valores xk del método implícito de Euler.
    xk = x0
    
    for i in t0:h:tf
        push!(listt,i) 
    end

    for j in 1:length(listt) #Se realiza un for donde le da un intervalo a j cuyo tamaño es el de 
     #la lista de t.
        push!(listx,xk) #Se guardan los valores xk en listx.
        g = - f(listt[j], xk)*h #Se define una función g que se usará en el método implícito.
        dg = 1 - f(listt[j], xk + h) + f(listt[j],xk) #Se escribe la derivada de la función g.
        xk = xk - g/dg #Se realiza la operación que determina los valores xk
    end
    
    return listt, listx #Se regresan la listas listt y listx.
end

export metodo_explicito_regla_punto_medio

"""documentación del método explícito de Euler (regla del punto medio)"""

function metodo_explicito_regla_punto_medio(f,x0,t0,tf,h) #Se crea una función que permite, 
    #mediante el método del punto medio, obtener la solución aproximada de una ecuación 
    #diferencial. Las entradas de esta función son las mismas que aquellas de la función 
    #para el método implícito de Euler.
    listt_medio = [] #Se crea un vector cuyo valor inicial es t0 y valor final tf donde la 
    #separación entre sus valores está dada por el valor de h.
        for i in t0:h/2:tf
            push!(listt_medio,i)
        end
    listx = [] #Se crea un vector en donde se guardan los valores xk del método de la regla 
    #del punto medio.
    xk = x0
        for j in 1:length(listt_medio) #Se realiza un for donde se le asigna un intervalo a j 
        #cuyo tamaño es la mitad del tamaño de la lista t.
        if iseven(j) 
            push!(listx,xk) #Se guarda el valor de j en la lista x.
            xk = xk + f(listt_medio[j+1],xk+ f(xk,listt_medio[j]*h/2))*h #Se crea una función 
            #que se requiere en este método.
        end
        
        end
    push!(listx,xk)
        listt = [] #Se crea un vector cuyo valor inicial es t0 y valor final tf donde la 
    #separación entre sus valores está dada por el valor de h.
        for i in t0:h:tf
            push!(listt,i)
        end
    return listt, listx #Se regresan la listas listt y listx.
end

export metodo_runge_kutta_vec

"""documentación del método de Runge-Kutta de orden 4 vectorial"""

function metodo_runge_kutta_vec(f,list,x0) #Se crea una función que permite, mediante el
    #método de Runge-Kutta de orden 4 vectorial, obtener la solución aproximada de una 
    #ecuación diferencial. Las entradas de esta función son la función, el valor inicial x0 
    #y una lista.
    
    x = x0
    h = list[2]-list[1] #Se calcula el tamaño del paso.
    listx = [] #Se crea un vector en donde se guardan los valores xk del método de Runge-Kutta de 
    #orden 4.
    push!(listx,x) #Se  guardan la los valores de x en listx.
    
        for i in 2:length(list) #Se realiza un for donde se da un intervalo a i cuyo tamaño es el
        #de la lista t.
            t = i*h
            #Se definen funciones necesarias para implementar el método de Runge-Kutta de orden 4.
            k1 = f(x,t)
            k2 = f(x + ((h/2) * k1),t + h/2)
            k3 = f(x + ((h/2) * k2),t + h/2)
            k4 = f(x + h * k3,t + h)
            x = x + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)#Se realiza la operación que determina los 
            #valores xk.
            push!(listx,x) #Se guarda el valor de xk en la lista x.
        end
     return listx #Se regresa la lista listx.
end

end
