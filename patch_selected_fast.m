

function [B,L]=patch_selected_fast(Aph,BTH,LBL,ths)

         ui=mean(Aph);
         uj=mean(BTH);

         term1= (2*ui*uj)./(ui^2+uj.^2);
         
         si=std(Aph).^2;
         ssj=std(BTH).^2;
         sisj=abs(1/(size(BTH,1)-1)*sum(repmat((Aph-ui),[1,size(BTH,2)]).*(BTH-repmat(uj,[size(BTH,1) 1]))));
         
         term2 = 2*sisj./(si+ssj);
         ss = term1.*term2;
         
         ssc = ss>ths;
         
        if sum(ssc)==0;
           ssc=(ss>(max(ss)*ths));
        end
        
        B=BTH(:,ssc);
        L=LBL(ssc);
         
end
