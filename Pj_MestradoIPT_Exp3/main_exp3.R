###########################################################################################
# Instituto de Pesquisas Tecnologicas de São Paulo                                        #
#                                                                                         #
# Aluno: Fernando de Lima Vicente                                                         #
# Curso: Infraestrutura Computacional                                                     #
#                                                                                         #
# Descrição: Testes com dados do dataset Intel Lab.                                       #
#            (Fonte 1: http://db.csail.mit.edu/labdata/data.txt.gz)                       #
#            (Fonte 2: http://db.csail.mit.edu/labdata/labdata.html)                      #
#                                                                                         #
#            Cabeçalho utilizado no dataset (Mica2Dot + Weather Board):                   #
#            date time epoch	moteid temperature humidity light voltage                   #
#                                                                                         #
#                                                                                         #
###########################################################################################

###########################################################################################
# Variáveis Globais                                                                       #
###########################################################################################

DATASET_FILE = "data.49.txt"

NODE = 49

QFKN = as.double(0)

RFKN = as.double(0)

                    # Limiar de temperatura (Qual a variação máxima que o experimento
TL = 5              # permite em uma iteração)

ZI = 2.9            # Controle de influência zero

GM_LEN = 0150       # Tamanho da janela do GM(1,1)

JI = 00000          # Janela inicial de observação

JF = 50000          # Janela final de observação

###########################################################################################
# Funções                                                                                 #
###########################################################################################


# Captura de dados ###########################################

obter_dados_dataset = function (node)                   {
  
  print ("Reading dataset file...")
  dataset = read.csv(DATASET_FILE, header = T, sep = " ", stringsAsFactors = F, dec=".")
  print ("Reading Complete....")
  
  print ("Filtering dataset by mote...")
  nodeset = dataset[which(dataset$moteid == node),]
  print ("Filtering Complete...")
  
  print ("Replacing NA by Zero...")
  nodeset$temperature [which(is.na(nodeset$temperature) == T)] = 0
  print ("Done...")
  
  return (nodeset)
  
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
  
  R = r            # Confiança no modelo matemático
  Q = 0.00075*q    # Confiança na leitura do sistema
  
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
    if( nodeVoltage[i] >= 2.8)
    {
      out[i] = 0.075
    }
    else
    {
      if(nodeVoltage[i] >= 2.7)
      {
        out[i] = 0.0075 
      }
      else
      {
        if(nodeVoltage[i] >= 2.6)
        {
          out[i] = 0.00075 
        }
        else
        {
          out[i] = 0.000075
        } 
      } 
    }  
  }
  
  return (out)
}

obter_confianca_node_old = function (nodeVoltage)       {
  
  out = double(0)
  
  for (i in 1:length(nodeVoltage))
  {
    if( nodeVoltage[i] >= 2.70)
    {
      out[i] = 0.001
    }
    else
    {
      if(nodeVoltage[i] >= 2.65)
      {
        out[i] = 0.0001 
      }
      else
      {
        if(nodeVoltage[i] >= 2.60)
        {
          out[i] = 0.0001 
        }
        else
        {
          out[i] = 0.0001
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
    
    R = RFKN[k] <<- 1                          # Confiança no modelo matemático
    
    Q = QFKN[k] <<- obter_confianca_node(c[k]) # Confiança na leitura do sistema
    
    
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
  
  for (i in 1:length (calculatedX1)) {
    calculatedX1 [i+1] = (X0[1] -1*(u/a))*exp (-1*a*i) + u/a
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
      # Precição                               #
      ##########################################
      
      out[win_start:win_end] = tail(predict_gm_model (c(tail(out,1), X0[win_start:win_end]),-1),win_size)
      
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

 #y = obter_dados_dataset(NODE)
 x = seq(0,nrow(y)-1)
 #x = x * 0.034

 r = q <- sd (y$temperature)

 
 #z = filtro_kalman_novo_3   ( y$temperature, y$voltage, r, q              )
 #w = filtro_kalman_classico ( y$temperature, r, q                         )
 #u = predict_gm_enhanced_2  ( y$temperature, length(y$temperature)/GM_LEN )
 
 #plot   ( x [JI:JF], y$voltage [JI:JF], type="l", col="black"  )
 plot   ( x [JI:JF], y$temperature [JI:JF], type="l", col="black", xlab="Tempo (Ks)",ylab="Temperatura (ºC)",ylim=c(0,40))
 #points ( x [JI:JF], w [JI:JF],             type="l", col="red"   )
 #points ( x [JI:JF], z [JI:JF],             type="l", col="red"    )
 #points ( x [JI:JF], u [JI:JF],             type="l", col="red" )

###########################################################################################
# Fim                                                                                     #
###########################################################################################
