using DelimitedFiles

include("GL.jl")

function main(num) #num指文件编号:num.txt, epoch_src_num.txt
	path="/Users/baotong/Desktop/period_Tuc/txt_all_obs_p90/"
    Tlist = readdlm(path*"$num.txt")[:,1] #这里、下面的epoch以及最后一行写入csv可以修改路径
    epoch = readdlm(path*"epoch_src_$num.txt")[:,1:2]
    gtis = [epoch[i,:] for i in 1:size(epoch,1)]
    ω_lo = 2*pi*1/50000 #频率下限
    ω_hi = 2*pi*1/10000 #频率上限
    dω = 2*pi*1e-8 #频率间隔
    r = compute_GL(Tlist, ω_lo:dω:ω_hi, gtis=gtis)
    writedlm(path*"result_$num.txt", [r[1] r[2] r[3] r[4]...], ',') #存入num.csv文件中，保存顺序与compute_GL一致
end

main(345)

#这里是编号从1到1000,如果一起做可以修改范围
#for i in 1:1000
#    main(i)
#end
