###########################################################################################
# Instituto de Pesquisas Tecnologicas de São Paulo                                        #
#                                                                                         #
# Aluno: Fernando de Lima Vicente                                                         #
# Curso: Infraestrutura Computacional                                                     #
#                                                                                         #
# Descrição: Testes com dados artificiais.                                                #
#                                                                                         #
###########################################################################################

###########################################################################################
# Variáveis Globais                                                                       #
###########################################################################################

PREC = 0.1           # Define a precisão da variavel "x"   

SIZE = 10000         # Define o tamanho da variável "x"

BATTERY_FULL = 3     # Tensão inicial da bateria

TL = 3               # Limiar de Temperatura 5

ZI = 2.9990          # Controle de influência zero

GM_LEN =900          # Tamanho da janela do GM(1,1)

QFKN = as.double(0)

RFKN = as.double(0)

JI = 000000          # Janela inicial de observação

JF = 100000          # Janela final de observação

###########################################################################################
# Funções                                                                                 #
###########################################################################################


# Geração de dados ###########################################

gerar_sensor_var = function ()                          {
  
  x = seq (0,SIZE, PREC)
  
  y = sin (x/2000) + 20
  
  return (y)
}

gerar_bateria_sensor = function ()                      {
  
  # Descarga da bateria
  
  x = seq (0,SIZE, PREC)
  y = numeric(0)
  
  # Função linear decrescente
  #y = -1 * (BATTERY_FULL / SIZE) * x + BATTERY_FULL
  
  # Função exponencial decrescente
  #y = BATTERY_FULL * (-1 * (11/10)^(x-SIZE) + 1)

  # Função exponencial decrescente
  y = BATTERY_FULL * (1 - ((1001/1000)^(x-SIZE)))
  
  y[1:2000] = BATTERY_FULL
  
  return (y)
}
  
gerar_vies_sensor_var = function (nodeVar, nodeVoltage) {
  
  vies = rnorm (SIZE / PREC + 1) 
  
  vp = nodeVoltage / BATTERY_FULL # Porcentagem da voltagem 
  
  out = nodeVar + (1 - vp) * vies # Começa em 0% até 100%
  
  return (out)
}


# Filtro de Kalman Clássico ##################################

filtro_kalman_classico = function (z, r, q)             {
  
  A = 1
  B = 0
  H = 1 
  I = 1
  
  x = z[1]
  u = 1:length(z)
  P = 1
  
  R = r           # Confiança no modelo matemático
  Q = 0.00000075*q    # Confiança na leitura do sistema
  
  xp = numeric(0)
  Pp = numeric(0)
  K  = numeric(0)
  
  print ("Processando Filtro de Kalman Classico: ")
  
  ##########################################
  # Status                                 #
  ##########################################
  total = length(z)
  last = 0
  print ("Status: 000%")
  
  for (k in 2:total){
    
    ##########################################
    # Status                                 #
    ##########################################
    status = as.integer(k/total*100)
    
    if (status %% 10 == 0 && status != last)
    {
      last = status
      print (sprintf("Status: %.3d%%", status))  
    }
    
    ##########################################
    # Prediction                             #
    ##########################################
    
    xp[k] = A*x[k-1] + B*u[k] 
    
    Pp[k] = A*P[k-1]*t(A) + Q
    
    ##########################################
    # Correction                             #
    ##########################################
    
    K[k] = (Pp[k]*t(H)) / (H*Pp[k]*t(H) + R)
    
    x[k] = xp[k] + K[k]*(z[k] - H*xp[k])
    
    P[k] = (I - K[k]*H)*Pp[k]
    
  }
  
  return (x);
}


# Filtro de Kalman Modificado ################################

obter_confianca_node = function (nodeVoltage)           {
  
  out = double(0)
  
  for (i in 1:length(nodeVoltage))
  {
    if( nodeVoltage[i] >= 2.99)
    {
      out[i] = 0.75
    }
    else
    {
      if(nodeVoltage[i] >= 2.9)
      {
        out[i] = 0.00075 #0.00001
      }
      else
      {
        if(nodeVoltage[i] >= 2.7)
        {
          out[i] = 0.0000075 #0.000001
        }
        else
        {
          out[i] = 0.00000075 #0.0000001
        } 
      } 
    }  
  }
  
  return (out)
}

filtro_kalman_novo = function (z, c)                    {
  
  A = 1
  B = 0
  H = 1 
  I = 1
  
  x = z[1]
  u = 1:length(z)
  P = 1
  
  R = RFKN[1] <<-  1                          # Confiança inicial no modelo matemático
  Q = QFKN[1] <<-  obter_confianca_node(c[1]) # Confiança inicial na leitura do sistema
  
  xp = numeric(0)
  Pp = numeric(0)
  K  = numeric(0)
  
  print ("Processando Filtro de Kalman Novo: ")
  
  ##########################################
  # Status                                 #
  ##########################################
  total = length(z)
  last = 0
  print ("Status: 000%")
  
  for (k in 2:total){
    
    ##########################################
    # Status                                 #
    ##########################################
    status = as.integer(k/total*100)
    
    if (status %% 10 == 0 && status != last)
    {
      last = status
      print (sprintf("Status: %.3d%%", status))  
    }
    
    ##########################################
    # Prediction                             #
    ##########################################
    
    xp[k] = A*x[k-1] + B*u[k] 
    
    Pp[k] = A*P[k-1]*t(A) + Q
    
    ##########################################
    # Correction                             #
    ##########################################
    
    K[k] = (Pp[k]*t(H)) / (H*Pp[k]*t(H) + R)
    
    x[k] = xp[k] + K[k]*(z[k] - H*xp[k])
    
    P[k] = (I - K[k]*H)*Pp[k]
  
    
    ##########################################
    # New Block                              #
    ##########################################
    
    R = RFKN[k] <<- 1.3                        # Confiança no modelo matemático
    
    Q = QFKN[k] <<- obter_confianca_node(c[k]) # Confiança na leitura do sistema
    
    
  }
  
  return (x);
}

filtro_kalman_novo_2 = function (z, c)                  {
  
  A = 1
  B = 0
  H = 1 
  I = 1
  
  x = z[1]
  u = 1:length(z)
  P = 1
  
  R = RFKN[1] <<-  sd(z)                              # Confiança inicial no modelo matemático
  Q = QFKN[1] <<-  sd(z) * obter_confianca_node(c[1]) # Confiança inicial na leitura do sistema
  
  xp = numeric(0)
  Pp = numeric(0)
  K  = numeric(0)
  
  #print ("Processando Filtro de Kalman Novo: ")
  
  ##########################################
  # Status                                 #
  ##########################################
  total = length(z)
  last = 0
  #print ("Status: 000%")
  
  for (k in 2:total){
    
    ##########################################
    # Status                                 #
    ##########################################
    status = as.integer(k/total*100)
    
    if (status %% 10 == 0 && status != last)
    {
      last = status
      #print (sprintf("Status: %.3d%%", status))  
    }
    
    ##########################################
    # Prediction                             #
    ##########################################
    
    xp[k] = A*x[k-1] + B*u[k] 
    
    Pp[k] = A*P[k-1]*t(A) + Q
    
    ##########################################
    # Correction (With outliar control)      #
    ##########################################
    
    K[k] = (Pp[k]*t(H)) / (H*Pp[k]*t(H) + R)
    
    ##########################################
    if (abs (z[k]) < (abs (x[k-1]) + TL)) {
      x[k] = xp[k] + K[k]*(z[k] - H*xp[k])  # Uses real value and estimation
    }
    else {
      x[k] = xp[k]                          # Uses estimation only
    }
    ##########################################
    
    P[k] = (I - K[k]*H)*Pp[k]
    
    
    ##########################################
    # New Block                              #
    ##########################################
    
    R = RFKN[k] <<- sd(z)                               # Confiança no modelo matemático
    
    Q = QFKN[k] <<- sd(z) * obter_confianca_node(c[k])  # Confiança na leitura do sistema
    
    
  }
  
  return (x);
}

filtro_kalman_novo_3 = function (z, c, r, q)            {
  
  A = 1
  B = 0
  H = 1 
  I = 1
  
  x = z[1]
  u = 1:length(z)
  P = 1
  
  R = RFKN[1] <<-  r                              # Confiança inicial no modelo matemático
  Q = QFKN[1] <<-  q * obter_confianca_node(c[1]) # Confiança inicial na leitura do sistema
  
  xp = numeric(0)
  Pp = numeric(0)
  K  = numeric(0)
  
  print ("Processando Filtro de Kalman Novo: ")
  
  ##########################################
  # Status Inicial                         #
  ##########################################
  total = length(z)
  last = 0
  print ("Status: 000%")
  
  for (k in 2:total){
    
    
    ##########################################
    # Status                                 #
    ##########################################
    status = as.integer(k/total*100)
    
    if (status %% 10 == 0 && status != last)
    {
      last = status
      print (sprintf("Status: %.3d%%", status))  
    }
    
    
    ##########################################
    # Predição                               #
    ##########################################
    
    xp[k] = A*x[k-1] + B*u[k] 
    
    Pp[k] = A*P[k-1]*t(A) + Q
    
    ##########################################
    # Correção com controle de outliers e    #
    # mecanismo de incluência zero           #
    ##########################################
    
    K[k] = (Pp[k]*t(H)) / (H*Pp[k]*t(H) + R)
    
    ##########################################
    if (c[k] >= ZI)
    {
      x[k] = z[k]                             # Não faz estimativa se a bateria está carregada
    }
    
    else
    {
      if (abs (z[k]) < (abs (x[k-1]) + TL)) {
        x[k] = xp[k] + K[k]*(z[k] - H*xp[k])  # Utiliza o valor estimado e o da leitura
      }
      else {
        x[k] = xp[k]                          # utiliza o valor estimado apenas 
      }
    }
    ##########################################
    
    P[k] = (I - K[k]*H)*Pp[k]
    
    
    ##########################################
    # Atualiza a confiança do sistema        #
    ##########################################
    
    R = RFKN[k] <<- r                               # Confiança no modelo matemático
    
    Q = QFKN[k] <<- q * obter_confianca_node(c[k])  # Confiança na leitura do sistema
    
    
    
  }
  
  return (x);
}


# Modelo Gray ou GM(1,1) #####################################

convert_x0_x1 = function (X0)                           {
  
  ##########################################
  # build AGO Sequence                     #
  ##########################################
  X1 = X0
  
  for (i in 1:(length(X1)-1)) {
    X1[i+1] = X1[i] + X1[i+1]
  }
  
  return (X1)
}

convert_x1_x0 = function (X1)                           {
  
  ##########################################
  # unbuild AGO Sequence                   #
  ##########################################
  X0 = X1
  
  for (i in 1:(length(X0)-1)) {
    X0[i+1] = X1[i+1] - X1[i]
  }
  
  return (X0)
}

get_b_matrix = function (X1)                            {
  
  B = matrix(nrow=length(X1)-1,ncol=2)
  
  for (i in 1:length(X1)-1){
    B[i,1]=-0.5*(X1[i]+X1[i+1])
    B[i,2]=1  
  }  
  
  #for (i in 1:length(X1)-1){
  #  B[i,1]=-0.5*(sum(X1[1:(i+1)]))
  #  B[i,2]=1  
  #}  
  
  return (B)
}

get_y_matrix = function (X0)                            {
  
  Y = matrix(X0[2:length(X0)], nrow=length(X0)-1, ncol=1)
  
  return (Y)
}


predict_gm_model = function (X0, predictions)           {
  
  ##########################################
  # build AGO Sequence                     #
  ##########################################
  
  X1 = convert_x0_x1 (X0)
  
  
  ##########################################
  # Solving dif eq model by least squares  #
  ##########################################
  
  B = get_b_matrix (X1)
  Y = get_y_matrix (X0)
  
  
  BtB = t(B) %*% B    # Transposed B Matrix plus B Matrix 
  iBtB = solve (BtB)  # Calcuate inversed Matrix from B
  BtY = t(B) %*% Y    # Transposed B Matrix plus Y Matrix
  
  
  # Parameters from Least Square Estimate  
  a = (iBtB %*% BtY) [1,1]
  u = (iBtB %*% BtY) [2,1]
  
  ##########################################
  # Perform prediction                     #
  ##########################################
  
  calculatedX1 = numeric (length(X1)+predictions)
  
  calculatedX1[1] = X1[1]
  
  if (a == 0 ){
    for (i in 1:length (calculatedX1)) {
      calculatedX1 [i+1] = X0[1] 
    }
    #return (calculatedX1)
  }
  else{
    for (i in 1:length (calculatedX1)) {
      calculatedX1 [i+1] = (X0[1] -1*(u/a))*exp (-1*a*i) + u/a
    }
    #return (convert_x1_x0 (calculatedX1))
  }
  
  return (convert_x1_x0 (calculatedX1))
  
}

predict_gm_enhanced = function (X0, win_size)           {
  
  win_size = as.integer(win_size)
  out = numeric(0)
  
  win_start = 1
  win_end = win_size
  
  while ((win_end + win_size) < length(X0)) {
    
    cat ("start=", win_start, " end=", win_end," out_len=", length(out),"/", length(X0), "\n")
    
    out[win_start:win_end] = predict_gm_model (X0[win_start:win_end],-1)         
    
    win_start = win_end + 1
    win_end   = win_end + win_size
  }
  
  cat ("startf=", win_start, " end=", win_end," out_len=", length(out), "\n")
  
  out[win_start:length(X0)] = predict_gm_model (X0[win_start:length(X0)],-1)         
  
  cat ("startf=", win_start, " end=", win_end," out_len=", length(out), "\n")
  
  
  
  return (out)
  
}

predict_gm_enhanced_2 = function (X0, win_size)         {
  
  win_size = as.integer (win_size)
  out = NULL
  
  win_start = 1
  win_end = win_size

  print ("Processando GM(1,1): ")
  
  
  ##########################################
  # Status Inicial                         #
  ##########################################
  total = length(X0)
  last = 0
  print ("Status: 000%")
  
  if (win_size > 0)
  {  
    while ((win_end + win_size) < length(X0)) {
      
      ##########################################
      # Status                                 #
      ##########################################
      status = as.integer(win_end/total*100)
      
      if (status %% 10 == 0 && status != last)
      {
        last = status
        print (sprintf("Status: %.3d%%", status))  
      }
      
      
      ##########################################
      # Precisão                               #
      ##########################################

      #print (sprintf("win_start: %d - win_end: %d", win_start, win_end))         
      out[win_start:win_end] = tail(predict_gm_model (c(tail(out,1), X0[win_start:win_end]),-1),win_size)

      #print(X0[win_start:win_end])
      #print(predict_gm_model(X0[win_start:win_end],-1))
      
      win_start = win_end + 1
      win_end   = win_end + win_size
    }
  }
  
  out[win_start:length(X0)] = predict_gm_model (X0[win_start:length(X0)],-1)         

  print ("Status: 100%")
    
  return (out)
  
}


###########################################################################################
# Rotina Principal                                                                        #
###########################################################################################

# Comente as três linhas abaixo para não gerar dados novos
#sensorVar     = gerar_sensor_var      ()
#sensorVoltage = gerar_bateria_sensor  ()  
#sensorViesVar = gerar_vies_sensor_var (sensorVar, sensorVoltage)

r = q <- sd (sensorViesVar) 

print ("Analise inicial: Comparando os 3 metodos no mesmo grafico...Pressione enter para iniciar")
#readline()

time = seq (0, SIZE, PREC)
 
fkn = filtro_kalman_novo_3(sensorViesVar, sensorVoltage, r, q)
fkc = filtro_kalman_classico(sensorViesVar, r, q)
gm2 = predict_gm_enhanced_2(sensorViesVar,as.integer(length(sensorViesVar)/GM_LEN))


plot   (time[JI:JF], sensorViesVar[JI:JF], type="l")
points (time[JI:JF], fkn[JI:JF], type="l", col="yellow")
points (time[JI:JF], fkc[JI:JF], type="l",col="blue")
points (time[JI:JF], gm2[JI:JF], type="l",col="green")
points (time[JI:JF], sensorVar[JI:JF], type="l", col="red")

#plot(time[JI:JF], sensorVoltage[JI:JF], type="l")

###########################################################################################
# Analise de eficiência de casa algoritmo (Quem está mais próximo do real)                #
###########################################################################################

print ("Analise de desempenho e histogramas... Pressione enter para continuar")
#readline()

absfkn = abs(sensorVar[JI:JF]-fkn[JI:JF])
absfkc = abs(sensorVar[JI:JF]-fkc[JI:JF])
absgm2 = abs(sensorVar[JI:JF]-gm2[JI:JF])

cc = 0 # contador do filtro de Kalman clássico
cn = 0 # contador do filtro de Kalman modificado 
cg = 0 # contador do Modelo Grey GM(1,1)

cc_t = 0
cn_t = 0
cg_t = 0

for (i in 1:length(absfkn)){

  absv = c(absfkn[i], absfkc[i], absgm2[i])
  
  sv = which(absv == min (absv))
  
  #print (sprintf("round %d - fkn= %f | fkc= %f | gm= %f | sv= %d ",i,absfkn[i], absfkc[i], absgm2[i], sv))
  
  if (1 %in% sv){
    cn = cn + 1
    
    if (absv[1]==absv[2]){
      cc = cc + 1
    }
    if (absv[1]==absv[3]){
      cg = cg + 1
    }
  } 
  
  if (2 %in% sv){
    cc = cc + 1
    
    if (absv[2]==absv[1]){
      cn = cn + 1
    }
    if (absv[2]==absv[3]){
      cg = cg + 1
    }
  } 

  if (3 %in% sv){
    cg = cg + 1
    
    if (absv[3]==absv[1]){
      cn = cn + 1
    }
    if (absv[3]==absv[2]){
      cc = cc + 1
    }
  }
  
  #########################################
  
  if ((i %% 10000) == 0){
    
    cn_t = cn_t + cn
    cg_t = cg_t + cg
    cc_t = cc_t + cc
    
    print(sprintf("i: %d |cn: %d | cc: %d | cg: %d ", i , cn, cc, cg))
    
    cn = 0
    cg = 0
    cc = 0
    
  }  

}

plot(time[JI:JF], absfkc, type="h")
plot(time[JI:JF], absfkn, type="h")
plot(time[JI:JF], absgm2, type="h")

print(sprintf("cn: %d | cc: %d | cg: %d " , cn_t, cc_t, cg_t))

###########################################################################################
# Fim                                                                                     #
###########################################################################################

