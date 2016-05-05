

##### Data Center Best_Policy_Algorithm
## Author: Ehsan Siavashi
## Date: Dec 15, 2015

# Please first download S.csv, S2.csv and S3.csv files into your directory. 
# See https://ehsansiavashi.shinyapps.io/ComNet-app for visiting the web application.

library(shiny)
library(shinyIncubator)
library(ggplot2)


server <- function(input, output) {
        output$inmatrix <- renderUI({
                matrixInput("Infl", "Select the detail Influence Matrix: (Aij)", 
                            as.data.frame(matrix(c(0.5,0.5,0,0.5,0.5,0,0.5,0.5,0),nrow=3,ncol=3, byrow = TRUE)))
        })
        output$mc1matrix <- renderUI({
                matrixInput("MC1", " MC matrix of computing nodes influencing (sending) no other nodes (Aii(0)): (row and column order: 'O'verload, 'N'ormal, 'U'nderload)", 
                            as.data.frame(matrix(c(0.5,0.4,0.1,0.25,0.5,0.25,0.15,0.4,0.45),nrow=3,ncol=3,byrow = TRUE)))
        })
        output$mc2matrix <- renderUI({
                matrixInput("MC2", " MC matrix of computing nodes influencing 1 node (Aii(1))", 
                            as.data.frame(matrix(c(0.4,0.45,0.15,0.2,0.45,0.35,0.12,0.33,0.55),nrow=3,ncol=3,byrow = TRUE)))
        })
        output$mc3matrix <- renderUI({
                matrixInput("MC3", " MC matrix of computing nodes influencing 2 nodes (Aii(2))", 
                            as.data.frame(matrix(c(0.3,0.5,0.2,0.15,0.4,0.45,0.1,0.25,0.65),nrow=3,ncol=3,byrow = TRUE)))
        })
        output$mc4matrix <- renderUI({
                matrixInput("MC4", " MC matrix of computing nodes influencing 3 nodes (Aii(3))", 
                            as.data.frame(matrix(c(0.2,0.55,0.25,0.1,0.35,0.55,0.03,0.22,0.75),nrow=3,ncol=3,byrow = TRUE)))
        })
        output$mc5matrix <- renderUI({
                matrixInput("MC5", " MC matrix of computing nodes influencing 4 nodes (Aii(4)", 
                            as.data.frame(matrix(c(0.1,0.60,0.30,0.05,0.30,0.65,0.01,0.14,0.85),nrow=3,ncol=3,byrow = TRUE)))
        })
        output$Smatrix <- renderUI({
                
                m<-c(rep(c(1,0,0),5),rep(c(0,1,0),5),rep(c(0,0,1),5),rep(c(1,0,0),5),rep(c(0,1,0),5),rep(c(0,0,1),5))
                
                matrixInput("S", " Matrix of the initial state of the network: (Col_order: 'O', 'N', 'U'; Row_orDer: Node-1...Node-n)", 
                            as.data.frame(matrix(m,nrow=30,ncol=3,byrow = TRUE)))
        })
        
        data<-eventReactive(input$action, {
                
                DCnum <-30
                
                Ss<-read.csv("www/S.csv")
                Ss<-as.matrix(Ss)
                Ss2<-read.csv("www/S2.csv")
                Ss2<-as.matrix(Ss2[1:30,2:30])
                Ss3<-read.csv("www/S3.csv")
                Ss3<-as.matrix(Ss3[1:30,2:31])

                
                if(input$select=="AveDeg: 3.63; Non-Symetric Adj-Mrtx"){
                        G<-as.matrix(Ss)  
                }else if(input$select=="AveDeg: 3.80; Semi-Symetric Adj-Mrtx"){
                        G<-as.matrix(Ss2)
                }else{
                        G<-as.matrix(Ss3)
                }
                
                D_t<- t(G) 
                D_logical<-as.logical(G)
                
                MC1<-input$MC1
                MC2<-input$MC2
                MC3<-input$MC3
                MC4<-input$MC4
                MC5<-input$MC5
                
                C_1 <- matrix(0, DCnum, DCnum)
                C_2 <- matrix(0, DCnum, DCnum)
                C_3 <- matrix(0, DCnum, DCnum)
                C_4 <- matrix(0, DCnum, DCnum)
                C_5 <- matrix(0, DCnum, DCnum)
                
                S0<-input$S 
                
                # the H_calculator function gets a condition matrix, a sending matrix and produces the H matrix for the next step.
                H_calculator <- function(input,M, S, C, D_t, sending, DCnum) {
                        E_t <- D_t * t(C) + diag(1,DCnum,DCnum) * (t(D_t) %*% (matrix(1,DCnum,DCnum)-t(C)))
                        
                        # Calculating A_0
                        A_0 <- cbind(M,matrix(0,3,3*DCnum-3))
                        
                        O <- matrix(0,3,3)
                        for (i in 2:DCnum) {
                                A<-O
                                for (j in 2:DCnum) {
                                        if(i==j && (sum(sending[i,])==1)){
                                                A<-cbind(A, E_t[i,i]*MC1)
                                        }
                                        else if(i==j && (sum(sending[i,])==2)){
                                                A<-cbind(A, E_t[i,i]*MC2)
                                        }
                                        else if(i==j && (sum(sending[i,])==3)){
                                                A<-cbind(A, E_t[i,i]*MC3)
                                        }
                                        else if(i==j && (sum(sending[i,])==4)){
                                                A<-cbind(A, E_t[i,i]*MC4)
                                        }
                                        else if(i==j && (sum(sending[i,])==5)){
                                                A<-cbind(A, E_t[i,i]*MC5)
                                        }
                                        else if(i==j && (sum(sending[i,])==6)){
                                                A<-cbind(A, E_t[i,i]*MC5)
                                        }
                                        else if(i==j && (sum(sending[i,])==7)){
                                                A<-cbind(A, E_t[i,i]*MC5)
                                        }
                                        else if(i==j && (sum(sending[i,])==8)){
                                                A<-cbind(A, E_t[i,i]*MC5)
                                        }
                                        else{
                                                A<-cbind(A,O)
                                        }
                                }
                                A_0 <- rbind(A_0,A)
                        }
                        H <- A_0 + kronecker(E_t-(diag(1, DCnum,DCnum)*E_t),input$Infl)
                        H
                }
                
                # the P_calculator function gets the current state and matrix H and produces the probability matrix.
                P_calculator <- function(S, H) {
                        P <- as.vector(t(S)) %*% H
                        P
                }
                
                # Realization function
                Realization <- function(DCnum, prob, S) {
                        
                        for (k in 0:(DCnum-1)) {
                                y<-runif(1)
                                if ((prob[3*k+1])>y){
                                        S[k+1,1]<-1
                                }
                                else if (((prob[3*k+1])+(prob[3*k+2]))<y){
                                        S[k+1,3]<-1
                                }
                                else{   S[k+1,2]<-1
                                }
                        }
                        S
                }
                
                NUM<-100
                
                ## SCENARIO 1
                Sum1=numeric()
                S<-S0
                for (variable in 1:NUM) {
                        # C1
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_1[i,j]<-1
                                        }else
                                                C_1[i,j] <- S[i,3]*S[j,1]
                                        
                                }
                        }
                        # sending ij is 1 iff i is sending workload to j, otherwise it is 0.
                        sending_1 <- t(C_1 * D_logical)
                        H1<-H_calculator(input,MC1, S, C_1, D_t, sending_1, DCnum)
                        p1<- as.vector(t(S)) %*% H1
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, p1, S)
                        
                        normalSum1<- p1[2]
                        
                        k<-5
                        while (k <=(90)) {
                                normalSum1 <- normalSum1 + p1[k]
                                k<-k+3
                        }
                        Sum1<-c(Sum1, normalSum1)
                }
                
                
                ## Scenario 2:
                Sum2<-numeric()
                S<-S0
                for (variable in 1:NUM) {
                        #C2
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_2[i,j]<-1
                                        }else
                                                C_2[i,j] <- S[i,3]*S[j,1] + S[i,3]*S[j,2]
                                }
                        }
                        sending_2 <- t(C_2 * D_logical)
                        H2<-H_calculator(input,MC1, S, C_2, D_t, sending_2, DCnum)
                        p2<- as.vector(t(S)) %*% H2
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, p2, S)
                        
                        normalSum2<- p2[2]
                        k<-5
                        while (k <= (90)) {
                                normalSum2 <- normalSum2 + p2[k]
                                k<-k+3
                        }
                        Sum2<-c(Sum2, normalSum2)
                        
                }
                
                
                # Scenario 3: 
                Sum3<-numeric()
                S<-S0
                for (variable in 1:NUM) {
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_3[i,j]<-1
                                        }else
                                                C_3[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1]
                                }
                        }
                        sending_3 <- t(C_3 * D_logical)
                        H3<-H_calculator(input,MC1, S, C_3, D_t, sending_3, DCnum)
                        p3<- as.vector(t(S)) %*% H3
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, p3, S)
                        
                        normalSum3<- p3[2]
                        k<-5
                        while (k <= (90)) {
                                normalSum3 <- normalSum3 + p3[k]
                                k<-k+3
                        }
                        Sum3<-c(Sum3, normalSum3)
                        
                }
                
                # Scenario 4: 
                Sum4<-numeric()
                S<-S0
                for (variable in 1:NUM) {
                        #C4
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_4[i,j]<-1
                                        }else
                                                C_4[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1] + S[i,3]*S[j,2]
                                }
                        }
                        sending_4 <- t(C_4 * D_logical)
                        H4<-H_calculator(input,MC1, S, C_4, D_t, sending_4, DCnum)
                        p4<- as.vector(t(S)) %*% H4
                        
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, p4, S)
                        
                        normalSum4<- p4[2]
                        k<-5
                        while (k <= (90)) {
                                normalSum4 <- normalSum4 + p4[k]
                                k<-k+3
                        }
                        Sum4<-c(Sum4, normalSum4)
                        
                }
                
                
                # Scenario 5: 
                Sum5<-numeric()
                S<-S0
                for (variable in 1:NUM) {
                        # C5
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_5[i,j]<-1
                                        }else
                                                C_5[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1] + S[i,3]*S[j,2] + S[i,2]*S[j,2]
                                }
                        }
                        
                        sending_5 <- t(C_5 * D_logical)
                        H5<-H_calculator(input,MC1, S, C_5, D_t, sending_5, DCnum)
                        p5<- as.vector(t(S)) %*% H5
                        
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, p5, S)
                        
                        normalSum5<- p5[2]
                        k<-5
                        while (k <= (90)) {
                                normalSum5 <- normalSum5 + p5[k]
                                k<-k+3
                        }
                        Sum5<-c(Sum5, normalSum5)
                }
                
                TABLE1 <- rbind(Sum1, Sum2, Sum3, Sum4,Sum5, rep(0,100))
                TABLE1/3000   #N-Expectancy
                
                
                #=============================================================
                
                Best_Sum<-0
                Sum<- matrix(c(0,0,0,0,0),nrow = 5,ncol = 1)
                
                
                S<-S0
                for (variable in 1:100) {
                        #C1
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_1[i,j]<-1
                                        }else
                                                C_1[i,j] <- S[i,3]*S[j,1]
                                        
                                }
                        }
                        #C2
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_2[i,j]<-1
                                        }else
                                                C_2[i,j] <- S[i,3]*S[j,1] + S[i,3]*S[j,2]
                                }
                        }
                        #C3
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_3[i,j]<-1
                                        }else
                                                C_3[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1]
                                }
                        }
                        #C4
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_4[i,j]<-1
                                        }else
                                                C_4[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1] + S[i,3]*S[j,2]
                                }
                        }
                        #C5
                        for (i in 1:DCnum) {
                                for (j in 1:DCnum) {
                                        if(i==j){
                                                C_5[i,j]<-1
                                        }else
                                                C_5[i,j] <- S[i,3]*S[j,1] + S[i,2]*S[j,1] + S[i,3]*S[j,2] + S[i,2]*S[j,2]
                                }
                        }
                        
                        # sending ij is 1 iff i is sending workload to j, otherwise it is 0.
                        
                        sending_1 <- t(C_1 * D_logical)
                        H1<-H_calculator(input,MC1, S, C_1, D_t, sending_1, DCnum)
                        p1<- as.vector(t(S)) %*% H1
                        
                        sending_2 <- t(C_2 * D_logical)
                        H2<-H_calculator(input,MC1, S, C_2, D_t, sending_2, DCnum)
                        p2<- as.vector(t(S)) %*% H2
                        
                        sending_3 <- t(C_3 * D_logical)
                        H3<-H_calculator(input,MC1, S, C_3, D_t, sending_3, DCnum)
                        p3<- as.vector(t(S)) %*% H3
                        
                        sending_4 <- t(C_4 * D_logical)
                        H4<-H_calculator(input,MC1, S, C_4, D_t, sending_4, DCnum)
                        p4<- as.vector(t(S)) %*% H4
                        
                        sending_5 <- t(C_5 * D_logical)
                        H5<-H_calculator(input,MC1, S, C_5, D_t, sending_5, DCnum)
                        p5<- as.vector(t(S)) %*% H5
                        
                        P_table<- rbind(p1,p2,p3,p4,p5)
                        
                        normalSum<- P_table[,2]
                        
                        k<-5
                        while (k < (90)) {
                                normalSum <- normalSum + P_table[,k]
                                k<-k+3
                        }
                        
                        Sum <- Sum + normalSum
                        
                        if(which.max(normalSum)==1){
                                P<-p1
                                Best_Sum <- Best_Sum + normalSum[1]
                        }else if(which.max(normalSum)==2){
                                P<-p2
                                Best_Sum <- Best_Sum + normalSum[2]
                        }else if(which.max(normalSum)==3){
                                P<-p3
                                Best_Sum <- Best_Sum + normalSum[3]
                        }else if(which.max(normalSum)==4){
                                P<-p4
                                Best_Sum <- Best_Sum + normalSum[4]
                        }else if(which.max(normalSum)==5){
                                P<-p5
                                Best_Sum <- Best_Sum + normalSum[5]
                        }
                        
                        
                        
                        # Realization
                        S<-matrix(0,DCnum, 3) 
                        S<- Realization(DCnum, P, S)
                        
                }
                
                TABLE2<-rbind(Sum,Best_Sum)/3000
                #=============================================================
                TABLE <- cbind(TABLE1,TABLE2)
                TABLE
        })
        
        output$table1 <- renderTable({
                dat<-data()
                dat<-as.data.frame(dat[1:6,101])
                row.names(dat) <- c("Policy 1","Policy 2","Policy 3","Policy 4","Policy 5","Best Policy")
                names(dat)<-c("N-Expectancy")
                dat
        })
        
        output$plot1<- renderPlot({
                d_f<-data()
                d_f<-as.data.frame(cbind(c("Policy 1","Policy 2","Policy 3","Policy 4","Policy 5","Best Policy"), d_f[1:6,101]))
                names(d_f)<-c("Policies", "N_Expectancy")
                
                ggplot(d_f, aes(x = factor(Policies), y = N_Expectancy)) + geom_bar(stat = "identity", fill="sky blue", colour = "black") 
        })
        
        output$table2 <- renderTable({
                n1=0
                n2=0
                n3=0
                n4=0
                n5=0
                for (i in 1:100) {
                        if(which.max(data()[,i])==1){
                                n1<-n1+1
                        }else if(which.max(data()[,i])==2){
                                n2<-n2+1
                        }else if(which.max(data()[,i])==3){
                                n3<-n3+1
                        }else if(which.max(data()[,i])==4){
                                n4<-n4+1
                        }else{
                                n5<-n5+1
                        }
                }
                dat<- as.matrix(c(n1,n2,n3,n4,n5))
                colnames(dat)<- c("Times Selected by the Best Policy")
                row.names(dat)<-c("Policy 1","Policy 2","Policy 3","Policy 4","Policy 5")
                dat
        })
        
        output$plot2<- renderPlot({
                n1=0
                n2=0
                n3=0
                n4=0
                n5=0
                for (i in 1:100) {
                        if(which.max(data()[,i])==1){
                                n1<-n1+1
                        }else if(which.max(data()[,i])==2){
                                n2<-n2+1
                        }else if(which.max(data()[,i])==3){
                                n3<-n3+1
                        }else if(which.max(data()[,i])==4){
                                n4<-n4+1
                        }else{
                                n5<-n5+1
                        }
                }
                T<-as.matrix(c(n1,n2,n3,n4,n5))
                T<-as.data.frame(T)
                T<-cbind(c("Policy 1","Policy 2", "Policy 3", "Policy 4", "Policy 5"),T)
                
                names(T)<-c("Policies","Times.Selected.by.Best.Policy")
                ggplot(T, aes(x = factor(Policies), y = Times.Selected.by.Best.Policy)) + geom_bar(stat = "identity", fill="sky blue", colour = "black") 

        })

}


ui <- fluidPage(
        
        theme = "bootstrap.min3.css",
        
        
        titlePanel(' Modeling a Network of Computing Nodes Using the Conditional Influence Model'),
        print(tags$em(HTML("<li> Created by Ehsan Siavashi</li>"))),
        
        tabsetPanel(
                tabPanel("Web application", 
                         
                         pageWithSidebar(
                                 print(" "),
                                 sidebarPanel(
                                         HTML("<p><b>If this is the first time You are using this app, 
                                              please read the documentations in the Documentation tab first.</b></p>
                                              <b>NOTE: Please do not change dimentions of the matrices.</b>" ),
                                         
                                         br(),

                                         selectInput("select", label = "Select a Directed Network of 30 Nodes (Influence Matrix D) with:", 
                                                     choices = c("AveDeg: 3.63; Non-Symetric Adj-Mtrx", 
                                                                 "AveDeg: 3.80; Semi-Symetric Adj-Mtrx",
                                                                 "AveDeg: 4.43; Symetric Adj-Mtrx" )),
                                         br(),
                                         uiOutput("inmatrix"),br(),
                                         
                                         uiOutput("mc1matrix"),
                                         uiOutput("mc2matrix"),
                                         uiOutput("mc3matrix"),
                                         uiOutput("mc4matrix"),
                                         uiOutput("mc5matrix"),
                                         uiOutput("Smatrix"),
                                         HTML("<br> &copy copyright researved."),
                                         wiDth = 4
                                 ),
                                 mainPanel(
                                         tags$hr(),
                                         print(tags$em(tags$strong("Set your desirable values in the sidebar panel, click the button below and wait for the simulation to complete (it may take a few seconds)."))),
                                         tags$hr(),
                                         actionButton(inputId = "action", label = "Run The Simulation"),
                                         tags$hr(),
                                         tableOutput("table1"),
                                         plotOutput("plot1"),
                                         tableOutput("table2"),
                                         plotOutput("plot2")
                                 )
                         )),
                tabPanel("Documentation",
                         
                         HTML("<p> <h3>The Conditional Influence Model</h3></p>
                              The conditional Influence Model (CIM) is an extension of the Influence Model (IM) introduced by C. Asavathiratham from MIT in 2001.
                              Influence Model is a tractable stochastic model for modeling dynamics of networked Markov Chains (MC). For more information about the influence model please visit <a href='https://dspace.mit.edu/handle/1721.1/33546'>DSpace@MIT</a>. </p>
                              <p>The Conditional Influence Model was introduced by Ehsan Siavashi and Dr. Mahshid Naeni from <a href='http://myweb.ttu.edu/marahnam/Files/ResearchGroup/DSNL.htm'>DSNL</a> lab in Texas Tech University and adds the following unique features: 
                                (1) The dynamics of components can vary in time, for instance, based on the internal on their past interactions, and (2) depending on the state of
                              the components the influences among components may get
                              activated or deactivated (rule-based interactions). For more information, link to the related publications will be available soon in this website. </p>
                              <p>Both IM and CIM are based on the idea that the next state of a node depends on its own MC and also the influences that each node receives from its neighbors.
                              </p>
                              <p><h3>Modeling Workload Balancing Policies in a Network of Computing Nodes Using CIM</h3></p>
                              <p> In this application we simulate several workload policies in a network of 30 computing nodes using CIM. The model also allows us to compare the different policies and suggest the best one. 
                                There are three possible states for each node: (O)verloaded, (N)ormal and (U)nderloaded.Node A influences Node B is interpreted as Node A sends workload to Node B.
                                The state transition probabilities are represented by a number of Markov chains. There are more than 1 MC for each node, because, depending on to how many neighbors a node is sending workload, 
                                its transition matrix differs. 
                                <p><h4> 1. Slider (Select the Network) </h4>
                              Allows the user to select among three different network topologies to examine the influence of these topologies on the results.
                                The networks are basically directed graphs which represent the Influence Matrix D introduced in the influence model (See <a href='https://dspace.mit.edu/handle/1721.1/33546'>here</a>.)
                              "),
                                br(),
                                img(height = 1200, width = 600, src="Final-plot.png"),
                                br(),
                                HTML("<p> <h4> 2. Markov Chain (MC) Matrices</h4>
                              <p>These four matrices represents the transition probabilities for each node based on the number of nodes it is sending workload to. Obviously a computing node that 
                                is not sending workload to any of its neighbors has a transition matrix that shows higher chance of going to a state with higher workload than a node which is sending its workload to a couple of its neighbors. 
                                Notice that in the default matrices the values in the columns from the left to the right columns increases from Aii(0) to Aii(4) (i.e. the chances of getting underloaded increases by sending workload to more neighbors). </p>
                              "),
                         br(),
                         img(height = 150, width = 275, src="z1.png"),
                         br(),
                         HTML("Example: The black square shows the chance of going to the O state assuming that the node is currently in the O state and
                              the red square shows the chance of going to the U state assuming that the node is currently in the N state. Notice that these matrices are row stochastic"),
                         
                         HTML("<h4> 3. Initial State (vector S[0])</h4>
                              <p>The matrix gets the initial state of the network. The matrix has 30 rows and three columns representing O, N and U states respectively. Therefore, the initial state of the nth node can be represented by setting 1 for its current state and 0 for the other two states. </p>"
                         ), 
                        HTML(("<h4> 4. Outputs</h4>
                               <p> (a) Table and Plot 1: These table and plot represent the N-expectancy of the network based on five different workload policies. N-expectancy is a measure for calculating the likelihood of going to N state in the next time step for a network. N-expectancy is a number between 0 and 1. For more information please refer to the paper. The policies are:<br>
                                <br>
                                Policy 1: Node i gets influenced by (receives workload
                                from) node j if and only if node i is U and
                              node j is O.<br><br>
                              Policy 2: node i gets influenced by node j if and
                              only if node i is U and node j is O or i is U and j
                              is N.<br><br>
                              Policy 3: node i gets influenced by node j if and
                              only if node i is U and node j is O or i is N and j
                              is O.<br><br>
                              Policy 4: node i gets influenced by node j if and
                              only if either node i is U and node j is O, i is U
                              and j is N or i is N and j is O.<br><br>
                                Policy 5: Node i gets influenced by node j if and only
                                if either node i is underloaded and node j is overloaded,
                              node i is underloaded and node j is in normal state, node
                              i is in normal state and node j is overloaded or node i
                              is in normal state and node j is in normal state too.<br><br>
                                The Best Policy is a policy that given the current state of the network selects the best policy amount Policy 1 to Policy 5 (the policy with highest N-expectancy). The results show that policy three is almost always the best policy amount the five policies above and also Policy 3 is almost as good as the Best policy.
                               <br>
                              <p> (b) Table and Plot 2: These table and matrix show the number of times that the Best Policy Algorithm selects each of the five policies in a run of 100 time steps. The results show that policy 3 is selected more than other policies and policy 2 (the worst policy) is very rarely selected.</p>")),
                         HTML("<br> &copy copyright researved. Ehsan Siavashi 2015")
                         )
                
                         )
        
                         )


shinyApp(ui = ui, server = server)

