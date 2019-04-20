library(ggtern);library(latex2exp)

rs = data.frame(pp = numeric(),p1 = numeric(),p2 = numeric(),a = numeric())

for (p_plas in seq(0.01,1,0.01)){
  p1 = seq(0.01,1-p_plas,0.01)
  p2 = 1-p1-p_plas
  pp = rep(p_plas,length(p1))
  rs <- rbind(rs,data.frame(pp,p1,p2))
}

finda0 <- function(r1,r2){
  hsc <- (1-r1-r2)^2
  a0  <- r1*hsc/(hsc*(1-r2))
  return(a0)
}

rs$a <- finda0(rs$p1,rs$p2)
rs[rs<0]=0
condition <- data.frame(p1=seq(0.35,0.55,by=0.01),p2=0.01,pp=0.59,a=0.5)

ggtern(rs,aes(p1,p2,pp,value=a))+
  geom_point(aes(color=a),alpha=0.05,stroke=0)+
  stat_interpolate_tern(aes(color=..level..),
                        breaks=seq(0,1,by=0.05))+
  scale_color_gradient2(low='darkblue',mid='green',high='darkred',
                        midpoint=mean(rs$a),
                        name=TeX('$$\\frac{u_{\\infty,2}}{u_{\\infty,1}}$$'))+
  theme_bw()+
  ggtitle(TeX('$$\\frac{u_{\\infty,2}}{u_{\\infty,1}}$$ needed for $$v_{abs,2}< 0$$'))+
  labs(x = TeX('$\\phi_{RBC}$'),
       y = TeX('$\\phi_{bac}$'),
       z = TeX('$\\phi_{plasma}$'))
  #geom_point(data=condition,aes(p1,p2,pp))

ggsave('C:/Users/cdhig/Documents/MS-Thesis-Pitt/Simulations/inversion_ternary.png')
