GB comments:
Prob1: 100%
Prob2:
P1:100
P2: 50 Did not find the longest ORF within frame of a start codon. 
P3:100
P4:100
P5:100 
Prob3
P1: 100
P2:100
P3:100 
Overall: 94


% Homework 1. Due before class on 9/5/17

%Sanjana Srinivasan

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
%x = 3; y = '5'; %mixed

%your code goes here
if isnumeric(x)
    x=x;
else
    x=str2num(x);
end

if isnumeric(y)
    y=y;
else
    y=str2num(y);
end

x+y
%output your answer
%x+y = 8
%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length
letters=["A", "T", "C", "G"];
nums = randi(numel(letters),[1 N]);
rand_seq = strjoin(letters(nums), "");


%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.
start=strfind(rand_seq, "ATG");
stop_codon1=strfind(rand_seq, "TAA");
stop_codon2=strfind(rand_seq, "TGA");
stop_codon3=strfind(rand_seq, "TAG");

stop_position=cat(2, stop_codon1, stop_codon2, stop_codon3)

all_orf=[]
for i = 1:length(start)
    orf=min(stop_position(stop_position > start(i)) - start(i));
    all_orf=[all_orf; orf];
end
max_orf=max(all_orf);
      
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
N = 500; % define sequence length
letters=["A", "T", "C", "G"];

count=0;
for i = 1:1000
    all_orf=[];
    nums = randi(numel(letters),[1 N]);
    rand_seq = strjoin(letters(nums), "");
    start=strfind(rand_seq, "ATG");
    stop_codon1=strfind(rand_seq, "TAA");
    stop_codon2=strfind(rand_seq, "TGA");
    stop_codon3=strfind(rand_seq, "TAG");
    stop_position=cat(2, stop_codon1, stop_codon2, stop_codon3);
    
    for j = 1:length(start)
        orf=min(stop_position(stop_position > start(j)) - start(j));
        all_orf=[all_orf; orf];
    end
    max_orf=max(all_orf);
    if max_orf >= 50
        count=count+1;
    end
end
prob_over50=count/1000;

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
mat=[]
N = randi([50,5000],1,100); % define sequence length
for v=1:length(N)
    letters=["A", "T", "C", "G"];
    count=0;
    for i = 1:1000
        all_orf=[];
        nums = randi(numel(letters),[1 N(v)]);
        rand_seq = strjoin(letters(nums), "");
        start=strfind(rand_seq, "ATG");
        stop_codon1=strfind(rand_seq, "TAA");
        stop_codon2=strfind(rand_seq, "TGA");
        stop_codon3=strfind(rand_seq, "TAG");
        stop_position=cat(2, stop_codon1, stop_codon2, stop_codon3);

        for j = 1:length(start)
            orf=min(stop_position(stop_position > start(j)) - start(j));
            all_orf=[all_orf; orf];
        end
        max_orf=max(all_orf);
        if max_orf >= 50
            count=count+1;
        end
    end
    prob_over50=count/1000
    mat1=[N(v) prob_over50]
    mat=[mat; mat1]
end

scatter(mat(:,1), mat(:,2));
title('Probability of Observing ORF > 50 bp');
xlabel('Length of Sequence');
ylabel('Probability');

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 
%SS: When N is very small, that is less than 50 base pairs, the probability 
%of having an ORF of 50 base pairs or more is 0. When N is higher than 50, 
%but still relatively small, if the DNA sequence is randomly generated as we have, 
%the probability of having a start and stop codon present would also be
%low, resulting in a low probability of having an ORF itself, and therefore
%a low probability of having one over 50 base pairs. As the length of the
%sequence increases, the probability of having at least one ORF increases, 
%and with increasing length of sequence and number of ORF's, the
%probability of obtaining ORF's longer than 50 base pairs also increases
%until it plateaus at a probability of 1.0. Thus, we would expect to see a
%hyperbolic curve with probability of 0.0 upto sequence length of 50, with
%probability increases with sequence length and plateauing at 1.0 at some
%length.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 
fid=fopen('qPCRdata.txt', 'r')
meta=fgetl(fid);
head=fgetl(fid);
line1=fgetl(fid);
line_vector=[];
while line1 ~= -1
    var1=strsplit(line1, "\t");
    var1str=cellfun(@str2num, var1, 'UniformOutput', false);
    line_vector=[line_vector var1str(5)];
    line1 = fgetl(fid);
end
fclose(fid);
line_vector = line_vector(:,1:end-24)

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
cp_matrix=cell2mat(line_vector)
cp_matrix=reshape(cp_matrix, [12,6]) 
cp_matrix=cp_matrix.'
% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

%taking the average expression for three replicates of each condition
for nn = 0:3    
    cp_avg(:,nn+1) = mean(cp_matrix(:,nn*3+1:(nn+1)*3),2); 
end


for i = (1:6) %rows
    for j = (1:3) %columns, excluding average expression of normalization gene in col 4
    norm_exp(i,j) = 2^(cp_avg(1,j) - cp_avg(i,j) - (cp_avg(1,4) - cp_avg(i,4)));
    end
end

bar(norm_exp)
title("2^ Normalized Gene Expression Across Six Conditions")
xlabel("Condition")
ylabel("2^ Normalized Cp")
legend("Gene1", "Gene2", "Gene3")

%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


