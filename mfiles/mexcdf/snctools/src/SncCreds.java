/* 
 * import org.apache.commons.httpclient.UsernamePasswordCredentials;
 */
import org.apache.commons.httpclient.*;
import org.apache.commons.httpclient.auth.*;

public class SncCreds implements CredentialsProvider {

    /* 
	 * Yeah, only username/password authentication for now.
	 */
    private UsernamePasswordCredentials myCredentials;

    public SncCreds(String userName, String password) {
        myCredentials = new UsernamePasswordCredentials(userName,password);
    }

	/*
	 * All we need do to implement the interface is implement
	 * the getCredentials method.
	 */
    public Credentials getCredentials (
		AuthScheme scheme, 
		String     host, 
		int        port, 
		boolean    proxy
	) throws CredentialsNotAvailableException {
		return(myCredentials);
	}

}
