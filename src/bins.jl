
#
# Main struct
#
struct NBin_data{T}

   # Bin ID
   ID::Int64

   # Average value
   μ::Vector{T}

   # weight do bin
   w::T

   # Variance
   s2::T

end

#
# Inputs
#
# X::Vector  -> a realization  ∈ R^nz
# Xm::Vector -> valores inferiores de todas as realizações ∈ R^nz
# Δ::Vector  -> "larguras" das divisões em cada dimensão ∈ R^nz
# Nb::Vector -> Número de divisões em cada dimensão ∈ R^nz
#
# Saida
#
# bin        -> número do bin que contém a realização
#
function Find_bin(X::Vector{T}, Xm::Vector{T}, Δ::Vector{T}, Nb::Vector{Int64}) where{T}

    # Verificação de consistência de dimensões
    length(X)==length(Xm)==length(Δ)==length(Nb) || throw("Find_Bin::dimensões devem ser iguais")
 
    # Verifica se as "larguras" são positivas
    all(Δ.>0) || throw("Find_Bin:: Δ deve ter todas as posições não nulas e positivas")

    # Primeira dimensão podemos fazer direto
    J = floor(Int64,(X[1]-Xm[1])/Δ[1]) 
    II = ifelse(J==Nb[1],J-1,J)

    # Primeira parte do cálculo do bin
    bin = II + 1
    
    # Loop pelas demais dimensões
    for i=2:length(X)
   
       # Indicador de posição nessa direção
       J = floor(Int64,(X[i]-Xm[i])/Δ[i]) 
       II = ifelse(J==Nb[i],J-1,J)

       # Acumula a posição do bin
       bin = bin + II*prod(Nb[1:i-1])
       
    end

    # Retorna o bin
    return bin

end


#
# Para cada realização em distrib, acha o bin associado e monta 
# a lista encadeada
#
# Entradas:
#
# distrib::Array -> matriz R^nz × R^n contendo as n realizações no R^nz
# Nb::Vector     -> Número de divisões em cada dimensão ∈ R^nz
#
# Saidas:
#
# next::Vector e head::Vector -> lista encadeada
#
function Linked_list(distrib::Array{T},Nb::Vector{Int64}) where {T}

   # Dimensões do problema
   nz, n = size(distrib)
   nbins = prod(Nb)
   
   # Aloca a lista encadeada
   head = zeros(Int64,nbins)  # posição em next
   next = -1*ones(Int64,n)    # número da realização

   # Valores mínimos, máximos e dimensões de cada bin
   xm = vec(minimum(distrib,dims=2))
   xM = vec(maximum(distrib,dims=2))
   Δ = (xM.-xm)./Nb
   
   # Loop sobre as realizações
   for i=1:n

      # realização
      X = vec(distrib[:,i])

      # Encontra o bin que contém a realização
      Ip = Find_bin(X, xm, Δ, Nb)

      # Encontra o ultimo ponteiro
      last_head = head[Ip]

      # Atualiza o ponteiro de bin
      head[Ip] = i

      # Se tem alguma realização já armazenada nesse bin, precisamos mudar o valor atual.
      # Mas antes, precisamos atualizar o ponteiro para a próxima posição
      if last_head > 0
         next[i] = last_head
      end   

   end #i

   # Retorna a lista encadeada
   return next, head

end

#
# Retorna a lista de realizações que estão em um dado bin
#
# Entradas:
#
# bin::Int64   -> um bin válido
# head e next  -> lista encadeada
#
# Saida:
#
# realizacoes_bin::Vector -> vetor contendo as realizações no bin
#
function Realizacoes_bin(bin::Int64, head::Vector{Int64}, next::Vector{Int64})

   # Verifica se bin é um valor válido
   0<bin<=length(head) || throw("Realizacoes_bin::bin $bin não é válido")
   
   # Vamos montar uma lista com as realizações que estão nesse bin
   realizacoes_bin = Int64[]

   # Acessa o ponteiro inicial para o bin
   pr = head[bin]

   # Vamos usar um sizehint para reduzir alocações
   hint = ceil(Int64,length(next)/10)
   sizehint!(realizacoes_bin, hint)

   while pr > 0

      # Guarda o número da realização
      push!(realizacoes_bin,pr)

      # Próxima posição na lista encadeada
      pr = next[pr]

   end

   return realizacoes_bin

end


#
# Gera uma lista com os dicionários de bins
#
function Generate_bins(dist::Array{T}, Nb::Vector{Int64};verbose=false) where {T}

   # Dimensão do problema
   nz, n = size(dist)

   !verbose || println("Generate_bins  nz=", nz, " n=",n)

   # Dicionário com os centros e pesos
   sampled = Dict{Int64,NBin_data}()

   # Verifica dimensões
   @assert nz==length(Nb) "Fatores deve ter dimensão == a $nz"
   @assert n>nz "Número de realizações deve ser maior do que a dimensão nz"

   # Vamos gerar a lista encadeada
   next,head = Linked_list(dist,Nb)

   # Loop pelos bins armazenando os bins efetivos
   contador = 1
   for i=1:prod(Nb)

      # Realizações para esse bin
      realizacoes = Realizacoes_bin(i, head, next)

      # Avalia se temos realizações nesse bin
      if !isempty(realizacoes)

         # Calcula média no bin
         μ = vec(mean(dist[:,realizacoes],dims=2))
         
         # Número de acessos
         na = length(realizacoes)

         # Peso
         w = na/n

         # Variação do bin 
         s2 = 0.0
         for r in realizacoes

            # Realização
            z = dist[:,r]

            # Realização menos a média
            dz = z .- μ

            # faz o produto interno e acumula
            s2 = s2 + dot(dz,dz)

         end #r   
         s2 = s2 / max((na-1),1)
         
         # Registra o bin
         sampled[contador] = NBin_data(i,μ,w,s2)
         contador += 1

      end

   end

   # Retorna o dicionário de bins efetivos
   return sampled

end 



