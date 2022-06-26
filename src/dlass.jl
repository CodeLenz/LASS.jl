#
# Calcula as derivadas da média e o desvio pelo LASS
#
# n é o número de amostras do MCS full 
#
function dLass(bins::Dict{Int64,NBin_data},f::Function,dfx::Function,nx::Int64,verbose = false)
   

    # Inicializa a derivada da média
    dmedia = zeros(nx)
    
    # Termos pra calcular a derivada da variância
    t1 = zeros(nx)
    
    # Vou aproveitar e calcular o valor esperado novamente aqui
    En = 0.0
    
    # Loop pelos bins para calcular a média e termos da variância
    for b in values(bins)
        
        # centro do bin
        xb = b.μ
        
        # Peso
        ω = b.w
       
        # Avalia a função no centro do bin e no valor da 
        # variável de projeto
        fxi = f(xb)
        
        # Avalia a derivada da função no centro do bin e no
        # valor da variável de projeto
        dfxi = dfx(xb)
        
        # Cálculo do valor esperado
        En += fxi*ω
        
        # Acumula a media da derivada
        dmedia .+= dfxi*ω
        
        # Acumula os termos para a variância
        t1 .+= fxi*dfxi*ω
        
              
    end #b
        
    # Calcula a derivada da variância
    dvariancia = 2*t1 .- 2*En*dmedia
    
    return dmedia, dvariancia

end
