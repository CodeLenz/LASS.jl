#
# Calcula a média e o desvio pelo MscMCS
#
function Lass(bins::Dict{Int64,NBin_data},f::Function, verbose = false)
   
    # Inicializa média
    media = 0.0
    
    # Termos pra calcular a variância
    t1 = 0.0
    
    # Loop pelos bins para calcular a média, já aproveitando
    # para calcular um dos termos da variância
    for b in values(bins)
        
        # centro do bin
        xb = b.μ
        
        # Peso
        ω = b.w
        
        # Avalia a função no centro do bin
        fxi = f(xb)
        
        # Acumula a media
        media += fxi*ω
        
        # Acumula o termo de E[f^2]
        t1 += ω*fxi^2
              
    end #b
    
    # Calcula o desvio
    variancia = t1 - media^2
    
    if verbose
       println("Média obtida pelo LASS ", media)
       println("Variância obtida pelo LASS ", variancia)
    end
    
    return media, variancia
end



